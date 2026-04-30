#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ete3 import Tree

from redip_utils import read_outgroup_species_from_genomes_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Root a gene tree on outgroup species specified in genomes.tsv."
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree.")
    parser.add_argument("--genomes-tsv", required=True, help="Pipeline genomes.tsv.")
    parser.add_argument("--output-tree", required=True, help="Output rooted Newick tree.")
    parser.add_argument("--summary-tsv", required=True, help="Rooting summary TSV.")
    parser.add_argument("--tip-separator", default="|", help="Tip-label separator.")
    parser.add_argument("--tree-format", type=int, default=1, help="ETE tree format.")
    return parser.parse_args()


def get_species(tip_name: str, separator: str) -> str:
    return tip_name.split(separator, 1)[0]


def root_on_present_outgroup(
    tree: Tree,
    outgroup_species: set[str],
    tip_separator: str,
) -> tuple[bool, str]:
    outgroup_leaves = [
        leaf for leaf in tree.iter_leaves()
        if get_species(leaf.name, tip_separator) in outgroup_species
    ]

    n_outgroup_tips = len(outgroup_leaves)

    if n_outgroup_tips == 0:
        return False, "no_outgroup_tips_present"

    found_species = sorted({
        get_species(leaf.name, tip_separator)
        for leaf in outgroup_leaves
    })

    if n_outgroup_tips == 1:
        tree.set_outgroup(outgroup_leaves[0])
        return True, "rooted_on_single_outgroup_tip__species=" + found_species[0]

    mrca = tree.get_common_ancestor(outgroup_leaves)

    if mrca is tree:
        return False, (
            "outgroup_mrca_is_root__tree_left_unchanged__n_outgroup_tips="
            + str(n_outgroup_tips)
            + "__species="
            + ",".join(found_species)
        )
    
    tree.set_outgroup(mrca)

    return True, (
        "rooted_on_outgroup_mrca__n_outgroup_tips="
        + str(n_outgroup_tips)
        + "__species="
        + ",".join(found_species)
    )


def write_summary(path: str, tree_id: str, status: str, message: str) -> None:
    clean_message = message.replace("\t", " ").replace("\n", " ")
    with open(path, "w") as out:
        out.write("tree_file\tstatus\tmessage\n")
        out.write(f"{tree_id}\t{status}\t{clean_message}\n")


def main() -> None:
    args = parse_args()

    tree_id = Path(args.tree).name
    outgroup_species = set(read_outgroup_species_from_genomes_tsv(args.genomes_tsv))

    try:
        tree = Tree(args.tree, format=args.tree_format)
    except Exception as exc:
        write_summary(
            args.summary_tsv,
            tree_id,
            "ERROR",
            f"failed_to_read_tree__{exc}",
        )
        sys.exit(1)

    try:
        rooted, message = root_on_present_outgroup(
            tree=tree,
            outgroup_species=outgroup_species,
            tip_separator=args.tip_separator,
        )

        tree.write(outfile=args.output_tree, format=args.tree_format)

        status = "ROOTED" if rooted else "SKIPPED"
        write_summary(args.summary_tsv, tree_id, status, message)

    except Exception as exc:
        write_summary(
            args.summary_tsv,
            tree_id,
            "ERROR",
            f"rooting_failed__{exc}",
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
