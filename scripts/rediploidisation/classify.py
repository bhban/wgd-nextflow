#!/usr/bin/env python3

"""
Adapted from classify_redip_history_of_focal_lineage_2026_02_12.py (Written by Drew Larson)
Classify rediploidisation histories in rooted gene trees.

This script supports three related use cases:

1. Standard mode
   The target species must have exactly --required-copies copies in a maximal
   allowed clade. This preserves the original classifier behaviour.

2. Recurrent-WGD balanced mode
   The target species may have more than --required-copies copies in a maximal
   allowed clade. The script reports complete recent duplicate groups and an
   ancestral event based on the MRCA of all target copies. If a clade has exactly
   --required-copies target copies, it is treated as standard_exact_n rather than
   recent_complete, because there is no evidence for a recurrent event.

3. Recurrent-WGD lossy mode
   If --allow-ancestral-loss is set, the script can report an ancestral event
   when there is at least one complete recent duplicate group plus at least one
   additional target-species copy, even if one expected recent group has been
   reduced by gene loss.

Example lossy topology for species tree (A,(B,C)):

    ((A,((B,C),(B,C))),(A,(B,C)));

For focal B and --required-copies 2, this can produce:

    recent     recent_complete          B1,B2
    recent     recent_singleton_lossy   B3
    ancestral  ancestral_lossy          B1,B2,B3

The ancestral row is classified using MRCA(B1,B2,B3).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import ete3

from redip_utils import read_redip_species


RecentGroup = Tuple[ete3.TreeNode, List[str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan rooted Newick trees, identify maximal allowed clades, and "
            "classify standard, recurrent, and lossy recurrent rediploidisation "
            "events for a focal target species."
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
        help=(
            "Copy-number rule used in standard mode. In recurrent mode, this "
            "is only used if --recurrent-wgd is not set."
        ),
    )
    parser.add_argument("--required-copies", type=int, default=2)

    parser.add_argument("--min-tips", type=int, default=1)
    parser.add_argument("--max-tips", type=int, default=None)

    parser.add_argument(
        "--recurrent-wgd",
        action="store_true",
        help=(
            "Allow more than --required-copies target copies in a maximal "
            "allowed clade. Reports recent complete duplicate groups and an "
            "ancestral event based on all target copies when criteria are met."
        ),
    )
    parser.add_argument(
        "--recent-grouping",
        choices=["auto", "monophyletic_exact", "closest_nonoverlapping"],
        default="auto",
        help=(
            "How to identify recent complete target-copy groups in recurrent "
            "mode. 'auto' tries monophyletic_exact first, then falls back to "
            "closest_nonoverlapping."
        ),
    )
    parser.add_argument(
        "--min-recent-groups",
        type=int,
        default=2,
        help=(
            "Minimum number of complete recent groups required to call a "
            "balanced ancestral recurrent event."
        ),
    )

    parser.add_argument(
        "--allow-ancestral-loss",
        action="store_true",
        help=(
            "In recurrent mode, allow ancestral events when one or more "
            "expected recent groups have been reduced by gene loss. Requires "
            "at least one complete recent group plus at least one additional "
            "target copy."
        ),
    )
    parser.add_argument(
        "--min-ancestral-target-copies",
        type=int,
        default=None,
        help=(
            "Minimum number of target-species copies needed to call a lossy "
            "ancestral recurrent event. Defaults to --required-copies + 1."
        ),
    )
    parser.add_argument(
        "--write-lossy-singletons",
        action="store_true",
        help=(
            "Write recent_singleton_lossy rows for unpaired target copies in "
            "lossy recurrent clades. These rows are useful for auditing but "
            "should not be used as Circos links."
        ),
    )

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


def parse_tip_label(label: str, separator: str, label_format: str) -> Dict[str, Optional[str]]:
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
    species_list: List[str] = []

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
    max_tips: Optional[int] = None,
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
    target_tips: List[str] = []

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


def target_tip_count(
    node: ete3.TreeNode,
    target_species: str,
    tip_separator: str,
    label_format: str,
) -> int:
    return len(
        find_target_tips(
            node=node,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )
    )


def find_monophyletic_exact_recent_groups(
    clade: ete3.TreeNode,
    target_species: str,
    tip_separator: str,
    label_format: str,
    group_size: int,
) -> List[RecentGroup]:
    groups: List[RecentGroup] = []

    for node in clade.traverse("postorder"):
        node_target_tips = find_target_tips(
            node=node,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        if len(node_target_tips) != group_size:
            continue

        child_has_same_target_count = False

        for child in node.children:
            child_count = target_tip_count(
                node=child,
                target_species=target_species,
                tip_separator=tip_separator,
                label_format=label_format,
            )

            if child_count == group_size:
                child_has_same_target_count = True
                break

        if child_has_same_target_count:
            continue

        groups.append((node, sorted(node_target_tips)))

    return remove_overlapping_groups(groups)


def remove_overlapping_groups(groups: Sequence[RecentGroup]) -> List[RecentGroup]:
    kept: List[RecentGroup] = []
    used: Set[str] = set()

    sorted_groups = sorted(groups, key=lambda x: (len(x[1]), ",".join(x[1])))

    for node, tips in sorted_groups:
        tip_set = set(tips)

        if tip_set & used:
            continue

        kept.append((node, tips))
        used.update(tip_set)

    return kept


def find_closest_nonoverlapping_recent_groups(
    full_tree: ete3.Tree,
    target_tips: Sequence[str],
    group_size: int,
) -> List[RecentGroup]:
    if group_size != 2:
        raise ValueError(
            "--recent-grouping closest_nonoverlapping currently supports "
            "--required-copies 2 only."
        )

    available = set(target_tips)
    pair_distances: List[Tuple[float, str, str]] = []

    for tip_a, tip_b in itertools.combinations(target_tips, 2):
        distance = full_tree.get_distance(tip_a, tip_b)
        pair_distances.append((distance, tip_a, tip_b))

    pair_distances.sort(key=lambda x: (x[0], x[1], x[2]))

    groups: List[RecentGroup] = []

    for _, tip_a, tip_b in pair_distances:
        if tip_a not in available or tip_b not in available:
            continue

        available.remove(tip_a)
        available.remove(tip_b)

        mrca = full_tree.get_common_ancestor([tip_a, tip_b])
        groups.append((mrca, sorted([tip_a, tip_b])))

    return groups


def find_recent_target_groups(
    full_tree: ete3.Tree,
    clade: ete3.TreeNode,
    target_species: str,
    tip_separator: str,
    label_format: str,
    group_size: int,
    recent_grouping: str,
) -> List[RecentGroup]:
    target_tips = find_target_tips(
        node=clade,
        target_species=target_species,
        tip_separator=tip_separator,
        label_format=label_format,
    )

    if len(target_tips) < group_size:
        return []

    if recent_grouping in {"auto", "monophyletic_exact"}:
        groups = find_monophyletic_exact_recent_groups(
            clade=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
            group_size=group_size,
        )

        if groups or recent_grouping == "monophyletic_exact":
            return groups

    if recent_grouping in {"auto", "closest_nonoverlapping"}:
        return find_closest_nonoverlapping_recent_groups(
            full_tree=full_tree,
            target_tips=target_tips,
            group_size=group_size,
        )

    raise ValueError(f"Unknown recent grouping mode: {recent_grouping}")


def append_event_row(
    rows: List[List[str]],
    tree_name: str,
    tree_index: int,
    event_level: str,
    event_type: str,
    loss_status: str,
    allowed_clade_index: int,
    num_tips_in_clade: int,
    num_target_tips_in_clade: int,
    event_index: int,
    target_tips: Sequence[str],
    sharing_species: Sequence[str],
) -> None:
    rows.append([
        tree_name,
        str(tree_index),
        event_level,
        event_type,
        loss_status,
        str(allowed_clade_index),
        str(num_tips_in_clade),
        str(num_target_tips_in_clade),
        str(event_index),
        ",".join(sorted(target_tips)),
        ",".join(sorted(sharing_species)),
    ])


def classify_tree_standard(
    tree: ete3.Tree,
    tree_name: str,
    tree_index: int,
    target_species: str,
    qualifying_clades: Sequence[ete3.TreeNode],
    tip_separator: str,
    label_format: str,
    copy_mode: str,
    required_copies: int,
) -> List[List[str]]:
    rows: List[List[str]] = []

    for clade_index, clade in enumerate(qualifying_clades, start=1):
        clade_labels = pull_tip_labels_from_ete3obj(clade)
        num_tips_in_clade = len(clade_labels)

        target_tips = find_target_tips(
            node=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        num_target_tips_in_clade = len(target_tips)

        if not clade_passes_copy_mode(
            node=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
            copy_mode=copy_mode,
            required_copies=required_copies,
        ):
            continue

        if len(target_tips) != required_copies:
            continue

        sharing_species = find_species_sharing_event(
            full_tree=tree,
            target_tips=target_tips,
            tip_separator=tip_separator,
            label_format=label_format,
        )

        append_event_row(
            rows=rows,
            tree_name=tree_name,
            tree_index=tree_index,
            event_level="standard",
            event_type="standard_exact_n",
            loss_status="complete",
            allowed_clade_index=clade_index,
            num_tips_in_clade=num_tips_in_clade,
            num_target_tips_in_clade=num_target_tips_in_clade,
            event_index=1,
            target_tips=target_tips,
            sharing_species=sharing_species,
        )

    return rows


def classify_tree_recurrent(
    tree: ete3.Tree,
    tree_name: str,
    tree_index: int,
    target_species: str,
    qualifying_clades: Sequence[ete3.TreeNode],
    tip_separator: str,
    label_format: str,
    required_copies: int,
    recent_grouping: str,
    min_recent_groups: int,
    allow_ancestral_loss: bool,
    min_ancestral_target_copies: int,
    write_lossy_singletons: bool,
) -> List[List[str]]:
    rows: List[List[str]] = []

    for clade_index, clade in enumerate(qualifying_clades, start=1):
        clade_labels = pull_tip_labels_from_ete3obj(clade)
        num_tips_in_clade = len(clade_labels)

        all_target_tips = sorted(
            find_target_tips(
                node=clade,
                target_species=target_species,
                tip_separator=tip_separator,
                label_format=label_format,
            )
        )

        num_target_tips_in_clade = len(all_target_tips)

        if num_target_tips_in_clade < required_copies:
            continue

        if num_target_tips_in_clade == required_copies:
            sharing_species = find_species_sharing_event(
                full_tree=tree,
                target_tips=all_target_tips,
                tip_separator=tip_separator,
                label_format=label_format,
            )

            append_event_row(
                rows=rows,
                tree_name=tree_name,
                tree_index=tree_index,
                event_level="standard",
                event_type="standard_exact_n",
                loss_status="complete",
                allowed_clade_index=clade_index,
                num_tips_in_clade=num_tips_in_clade,
                num_target_tips_in_clade=num_target_tips_in_clade,
                event_index=1,
                target_tips=all_target_tips,
                sharing_species=sharing_species,
            )

            continue

        recent_groups = find_recent_target_groups(
            full_tree=tree,
            clade=clade,
            target_species=target_species,
            tip_separator=tip_separator,
            label_format=label_format,
            group_size=required_copies,
            recent_grouping=recent_grouping,
        )

        used_recent_tips: Set[str] = set()

        for event_index, (_, group_target_tips) in enumerate(recent_groups, start=1):
            used_recent_tips.update(group_target_tips)

            sharing_species = find_species_sharing_event(
                full_tree=tree,
                target_tips=group_target_tips,
                tip_separator=tip_separator,
                label_format=label_format,
            )

            append_event_row(
                rows=rows,
                tree_name=tree_name,
                tree_index=tree_index,
                event_level="recent",
                event_type="recent_complete",
                loss_status="complete",
                allowed_clade_index=clade_index,
                num_tips_in_clade=num_tips_in_clade,
                num_target_tips_in_clade=num_target_tips_in_clade,
                event_index=event_index,
                target_tips=group_target_tips,
                sharing_species=sharing_species,
            )

        unused_target_tips = sorted(set(all_target_tips) - used_recent_tips)

        has_complete_recent_group = len(recent_groups) >= 1
        has_additional_focal_copy = len(unused_target_tips) >= 1

        exact_balanced_copy_count = (
            len(recent_groups) >= min_recent_groups
            and num_target_tips_in_clade == len(recent_groups) * required_copies
        )

        has_balanced_ancestral_signal = (
            len(recent_groups) >= min_recent_groups
            and num_target_tips_in_clade > required_copies
            and exact_balanced_copy_count
        )

        has_lossy_ancestral_signal = (
            allow_ancestral_loss
            and has_complete_recent_group
            and has_additional_focal_copy
            and num_target_tips_in_clade >= min_ancestral_target_copies
        )

        should_write_ancestral = (
            has_balanced_ancestral_signal
            or has_lossy_ancestral_signal
        )

        if write_lossy_singletons and has_lossy_ancestral_signal:
            singleton_start_index = len(recent_groups) + 1

            for offset, singleton_tip in enumerate(unused_target_tips, start=0):
                append_event_row(
                    rows=rows,
                    tree_name=tree_name,
                    tree_index=tree_index,
                    event_level="recent",
                    event_type="recent_singleton_lossy",
                    loss_status="lossy",
                    allowed_clade_index=clade_index,
                    num_tips_in_clade=num_tips_in_clade,
                    num_target_tips_in_clade=num_target_tips_in_clade,
                    event_index=singleton_start_index + offset,
                    target_tips=[singleton_tip],
                    sharing_species=[target_species],
                )

        if should_write_ancestral:
            if has_balanced_ancestral_signal and not has_lossy_ancestral_signal:
                ancestral_event_type = "ancestral_balanced"
                ancestral_loss_status = "complete"
            elif has_balanced_ancestral_signal and has_lossy_ancestral_signal:
                ancestral_event_type = "ancestral_lossy"
                ancestral_loss_status = "lossy"
            else:
                ancestral_event_type = "ancestral_lossy"
                ancestral_loss_status = "lossy"

            sharing_species = find_species_sharing_event(
                full_tree=tree,
                target_tips=all_target_tips,
                tip_separator=tip_separator,
                label_format=label_format,
            )

            append_event_row(
                rows=rows,
                tree_name=tree_name,
                tree_index=tree_index,
                event_level="ancestral",
                event_type=ancestral_event_type,
                loss_status=ancestral_loss_status,
                allowed_clade_index=clade_index,
                num_tips_in_clade=num_tips_in_clade,
                num_target_tips_in_clade=num_target_tips_in_clade,
                event_index=1,
                target_tips=all_target_tips,
                sharing_species=sharing_species,
            )

    return rows


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
    max_tips: Optional[int],
    recurrent_wgd: bool,
    recent_grouping: str,
    min_recent_groups: int,
    allow_ancestral_loss: bool,
    min_ancestral_target_copies: int,
    write_lossy_singletons: bool,
) -> List[List[str]]:
    qualifying_clades = find_maximal_allowed_clades(
        tree=tree,
        allowed_species=allowed_species,
        tip_separator=tip_separator,
        label_format=label_format,
        min_tips=min_tips,
        max_tips=max_tips,
    )

    if recurrent_wgd:
        return classify_tree_recurrent(
            tree=tree,
            tree_name=tree_name,
            tree_index=tree_index,
            target_species=target_species,
            qualifying_clades=qualifying_clades,
            tip_separator=tip_separator,
            label_format=label_format,
            required_copies=required_copies,
            recent_grouping=recent_grouping,
            min_recent_groups=min_recent_groups,
            allow_ancestral_loss=allow_ancestral_loss,
            min_ancestral_target_copies=min_ancestral_target_copies,
            write_lossy_singletons=write_lossy_singletons,
        )

    return classify_tree_standard(
        tree=tree,
        tree_name=tree_name,
        tree_index=tree_index,
        target_species=target_species,
        qualifying_clades=qualifying_clades,
        tip_separator=tip_separator,
        label_format=label_format,
        copy_mode=copy_mode,
        required_copies=required_copies,
    )


def prepare_output_file(path: Path, on_exists: str) -> str:
    if path.exists():
        if on_exists == "error":
            raise FileExistsError(f"Output file already exists: {path}")

        if on_exists == "overwrite":
            return "w"

        if on_exists == "append":
            return "a"

    return "w"


def output_header() -> List[str]:
    return [
        "tree_file",
        "tree_index",
        "event_level",
        "event_type",
        "loss_status",
        "allowed_clade_index",
        "num_tips_in_allowed_clade",
        "num_target_tips_in_allowed_clade",
        "event_index",
        "target_tips",
        "sharing_species",
    ]


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
            writer.writerow(output_header())

        for row in rows:
            writer.writerow(row)

    return len(rows)


def resolve_min_ancestral_target_copies(args: argparse.Namespace) -> int:
    if args.min_ancestral_target_copies is not None:
        return args.min_ancestral_target_copies

    return args.required_copies + 1


def validate_args(args: argparse.Namespace) -> None:
    if args.required_copies < 1:
        raise ValueError("--required-copies must be at least 1.")

    if args.min_recent_groups < 1:
        raise ValueError("--min-recent-groups must be at least 1.")

    min_ancestral_target_copies = resolve_min_ancestral_target_copies(args)

    if min_ancestral_target_copies <= args.required_copies:
        raise ValueError(
            "--min-ancestral-target-copies must be greater than "
            "--required-copies."
        )

    if args.allow_ancestral_loss and not args.recurrent_wgd:
        raise ValueError("--allow-ancestral-loss requires --recurrent-wgd.")

    if args.write_lossy_singletons and not args.allow_ancestral_loss:
        raise ValueError("--write-lossy-singletons requires --allow-ancestral-loss.")


def process_treefile(args: argparse.Namespace) -> List[List[str]]:
    all_rows: List[List[str]] = []
    tree_path = Path(args.treefile)
    min_ancestral_target_copies = resolve_min_ancestral_target_copies(args)

    allowed_species = set(
        read_redip_species(
            genomes_tsv=args.genomes_tsv,
            allowed_species=args.allowed_species,
        )
    )

    with open(args.treefile, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            if "\t" in line:
                tree_name, newick = line.split("\t", 1)
            else:
                tree_name = tree_path.name
                newick = line

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
                    recurrent_wgd=args.recurrent_wgd,
                    recent_grouping=args.recent_grouping,
                    min_recent_groups=args.min_recent_groups,
                    allow_ancestral_loss=args.allow_ancestral_loss,
                    min_ancestral_target_copies=min_ancestral_target_copies,
                    write_lossy_singletons=args.write_lossy_singletons,
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
    validate_args(args)

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

