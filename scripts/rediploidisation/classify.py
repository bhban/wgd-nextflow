#!/usr/bin/env python3

'''
Adapted from classify_redip_history_of_focal_lineage_2026_02_12.py (Written by Drew Larson)

Find maximal clades in gene trees whose species are all drawn from a user-supplied
allowed-species set, then report cases where the target species has the required
copy pattern and list the species descending from the MRCA of those copies.
'''

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set

import ete3

from redip_utils import read_redip_species


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan one or more Newick trees, identify maximal clades whose species "
            "are all drawn from the rediploidisation/WGD ingroup, and report "
            "target-species copy groups plus species sharing their inferred event."
        )
    )

    parser.add_argument("--target-species", required=True)
    parser.add_argument("--treefile", required=True)
    parser.add_argument("--genomes-tsv", default=None)
    parser.add_argument("--allowed-species", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tip-separator", default="@")
    parser.add_argument(
        "--label-format",
        choices=["species_gene", "species_chr_gene"],
        default="species_gene",
    )
    parser.add_argument(
        "--copy-mode",
        choices=[
            "target_exactly_n",
            "all_exactly_n",
            "target_exactly_n_others_min_n",
        ],
        default="target_exactly_n",
    )
    parser.add_argument("--required-copies", type=int, default=2)
    parser.add_argument("--min-tips", type=int, default=1)
    parser.add_argument("--max-tips", type=int, default=None)
    parser.add_argument(
        "--on-exists",
        choices=["overwrite", "append", "error"],
        default="overwrite",
    )
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--no-header", action="store_true")
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )

    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s: %(message)s",
    )


def parse_tip_label(label: str, separator: str, label_format: str) -> Dict[str, str | None]:
    parts = label.split(separator)

    if label_format == "species_gene":
        if len(parts) != 2:
            raise ValueError(
                f"Tip label '{label}' does not match format 'species_gene' "
                f"with separator '{separator}'."
            )

        species, gene_id = parts

        if not species or not gene_id:
            raise ValueError(f"Tip label '{label}' contains an empty field.")

        return {"species": species, "chr": None, "gene": gene_id}

    if label_format == "species_chr_gene":
        if len(parts) != 3:
            raise ValueError(
                f"Tip label '{label}' does not match format 'species_chr_gene' "
                f"with separator '{separator}'."
            )

        species, chr_id, gene_id = parts

        if not species or not chr_id or not gene_id:
            raise ValueError(f"Tip label '{label}' contains an empty field.")

        return {"species": species, "chr": chr_id, "gene": gene_id}

    raise ValueError(f"Unknown label format: {label_format}")


def pull_tip_labels_from_ete3obj(tree_node: ete3.TreeNode) -> List[str]:
    leaves: List[str] = []
    seen = set()

    for leaf in tree_node.iter_leaves():
        if leaf.name in seen:
            raise ValueError(f"Duplicate tip label detected: '{leaf.name}'")

        seen.add(leaf.name)
        leaves.append(leaf.name)

    return leaves


def get_species_from_labels(
    labels: Sequence[str],
    separator: str,
    label_format: str,
) -> List[str]:
    species_list = []

    for label in labels:
        parsed = parse_tip_label(label, separator, label_format)
        species = parsed["species"]

        if species is None:
            raise ValueError(f"Could not extract species from tip label '{label}'.")

        species_list.append(species)

    return species_list


def handle_error(message: str, strict: bool) -> None:
    if strict:
        raise ValueError(message)

    logging.warning(message)


def is_allowed_clade(
    node: ete3.TreeNode,
    allowed_species: Set[str],
    tip_separator: str,
    label_format: str,
) -> bool:
    labels = pull_tip_labels_from_ete3obj(node)
    species_set = set(get_species_from_labels(labels, tip_separator, label_format))
    return species_set.issubset(allowed_species)


def get_species_counts(
    node: ete3.TreeNode,
    tip_separator: str,
    label_format: str,
) -> Counter:
    labels = pull_tip_labels_from_ete3obj(node)
    species_list = get_species_from_labels(labels, tip_separator, label_format)
    return Counter(species_list)


def clade_passes_copy_mode(
    node: ete3.TreeNode,
    target_species: str,
    tip_separator: str,
    label_format: str,
    copy_mode: str,
    required_copies: int,
) -> bool:
    counts = get_species_counts(node, tip_separator, label_format)
    target_count = counts.get(target_species, 0)

    if required_copies < 1:
        raise ValueError("--required-copies must be at least 1.")

    if copy_mode == "target_exactly_n":
        return target_count == required_copies

    if copy_mode == "all_exactly_n":
        return bool(counts) and all(
            count == required_copies
            for count in counts.values()
        )

    if copy_mode == "target_exactly_n_others_min_n":
        if target_count != required_copies:
            return False

        for species, count in counts.items():
            if species == target_species:
                continue

            if count < required_copies:
                return False

        return True

    raise ValueError(f"Unknown copy mode: {copy_mode}")


def find_maximal_allowed_clades(
    tree: ete3.Tree,
    allowed_species: Set[str],
    tip_separator: str,
    label_format: str,
    min_tips: int = 1,
    max_tips: int | None = None,
) -> List[ete3.TreeNode]:
    maximal_nodes: List[ete3.TreeNode] = []

    for node in tree.traverse("postorder"):
        labels = pull_tip_labels_from_ete3obj(node)
        n_tips = len(labels)

        if n_tips < min_tips:
            continue

        if max_tips is not None and n_tips > max_tips:
            continue

        current_ok = is_allowed_clade(
            node=node,
            allowed_species=allowed_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        if not current_ok:
            continue

        parent = node.up

        if parent is None:
            maximal_nodes.append(node)
            continue

        parent_ok = is_allowed_clade(
            node=parent,
            allowed_species=allowed_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        if not parent_ok:
            maximal_nodes.append(node)

    return maximal_nodes


def find_target_tips(
    node: ete3.TreeNode,
    target_species: str,
    tip_separator: str,
    label_format: str,
) -> List[str]:
    labels = pull_tip_labels_from_ete3obj(node)
    target_tips = []

    for label in labels:
        parsed = parse_tip_label(label, tip_separator, label_format)

        if parsed["species"] == target_species:
            target_tips.append(label)

    return target_tips


def find_species_sharing_event(
    full_tree: ete3.Tree,
    target_tips: Sequence[str],
    tip_separator: str,
    label_format: str,
) -> List[str]:
    if len(target_tips) < 2:
        raise ValueError("At least two target tips are required to compute an MRCA.")

    mrca = full_tree.get_common_ancestor(list(target_tips))
    mrca_labels = pull_tip_labels_from_ete3obj(mrca)

    return sorted(set(get_species_from_labels(mrca_labels, tip_separator, label_format)))


def classify_tree(
    tree: ete3.Tree,
    tree_name: str,
    tree_index: int,
    target_species: str,
    allowed_species: Set[str],
    tip_separator: str,
    label_format: str,
    copy_mode: str,
    required_copies: int,
    min_tips: int,
    max_tips: int | None,
) -> List[List[str]]:
    rows: List[List[str]] = []

    qualifying_clades = find_maximal_allowed_clades(
        tree=tree,
        allowed_species=allowed_species,
        tip_separator=tip_separator,
        label_format=label_format,
        min_tips=min_tips,
        max_tips=max_tips,
    )

    for clade_index, clade in enumerate(qualifying_clades, start=1):
        clade_labels = pull_tip_labels_from_ete3obj(clade)
        num_tips_in_clade = len(clade_labels)

        if not clade_passes_copy_mode(
            node=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
            copy_mode=copy_mode,
            required_copies=required_copies,
        ):
            continue

        target_tips = find_target_tips(
            node=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        if len(target_tips) != required_copies:
            continue

        sharing_species = find_species_sharing_event(
            full_tree=tree,
            target_tips=target_tips,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        rows.append([
            tree_name,
            str(tree_index),
            str(clade_index),
            str(num_tips_in_clade),
            ",".join(target_tips),
            ",".join(sharing_species),
        ])

    return rows


def prepare_output_file(path: Path, on_exists: str) -> str:
    if path.exists():
        if on_exists == "error":
            raise FileExistsError(f"Output file already exists: {path}")

        if on_exists == "overwrite":
            return "w"

        if on_exists == "append":
            return "a"

    return "w"


def write_rows(
    output_path: str,
    rows: Iterable[List[str]],
    write_header: bool,
    mode: str,
) -> int:
    rows = list(rows)

    should_write_header = False
    if write_header:
        if mode == "w":
            should_write_header = True
        elif mode == "a":
            output_file = Path(output_path)
            if (not output_file.exists()) or output_file.stat().st_size == 0:
                should_write_header = True

    with open(output_path, mode, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        if should_write_header:
            writer.writerow([
                "tree_file",
                "tree_index",
                "allowed_clade_index",
                "num_tips_in_allowed_clade",
                "target_tips",
                "sharing_species",
            ])

        for row in rows:
            writer.writerow(row)

    return len(rows)


def process_treefile(args: argparse.Namespace) -> List[List[str]]:
    all_rows: List[List[str]] = []
    tree_path = Path(args.treefile)
    tree_name = tree_path.name

    allowed_species = set(
        read_redip_species(
            genomes_tsv=args.genomes_tsv,
            allowed_species=args.allowed_species,
        )
    )

    with open(args.treefile, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            newick = line.strip()

            if not newick:
                continue

            try:
                tree = ete3.Tree(newick)
            except Exception as exc:
                handle_error(
                    f"Could not parse Newick on line {line_number} of "
                    f"{args.treefile}: {exc}",
                    args.strict,
                )
                continue

            try:
                rows = classify_tree(
                    tree=tree,
                    tree_name=tree_name,
                    tree_index=line_number,
                    target_species=args.target_species,
                    allowed_species=allowed_species,
                    tip_separator=args.tip_separator,
                    label_format=args.label_format,
                    copy_mode=args.copy_mode,
                    required_copies=args.required_copies,
                    min_tips=args.min_tips,
                    max_tips=args.max_tips,
                )
                all_rows.extend(rows)

            except Exception as exc:
                handle_error(
                    f"Failed while processing tree on line {line_number} of "
                    f"{args.treefile}: {exc}",
                    args.strict,
                )
                continue

    return all_rows


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    output_path = Path(args.output)
    mode = prepare_output_file(output_path, args.on_exists)

    rows = process_treefile(args)

    write_rows(
        output_path=args.output,
        rows=rows,
        write_header=not args.no_header,
        mode=mode,
    )

    logging.info("Finished.")
    logging.info("Output file: %s", args.output)
    logging.info("Rows written: %d", len(rows))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.error(str(exc))
        sys.exit(1)
