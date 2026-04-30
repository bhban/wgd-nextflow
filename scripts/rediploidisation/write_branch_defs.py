#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Set

from ete3 import Tree

from redip_utils import read_redip_species


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate branch-definition TSV files for each focal species in the "
            "rediploidisation/WGD ingroup."
        )
    )
    parser.add_argument("--tree", required=True)
    parser.add_argument("--genomes-tsv", default=None)
    parser.add_argument("--allowed-species", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tree-format", type=int, default=1)
    return parser.parse_args()


def read_tree(path: str, tree_format: int) -> Tree:
    tree_text = Path(path).read_text().strip()

    if not tree_text:
        raise ValueError(f"Tree file is empty: {path}")

    return Tree(tree_text, format=tree_format)


def get_tree_tip_names(tree: Tree) -> Set[str]:
    return {leaf.name for leaf in tree.iter_leaves()}


def get_allowed_descendants(node, allowed_species_set: Set[str]) -> Set[str]:
    return {
        leaf.name
        for leaf in node.iter_leaves()
        if leaf.name in allowed_species_set
    }


def make_branch_sets_for_focal(
    tree: Tree,
    focal_species: str,
    allowed_species: List[str],
) -> List[Set[str]]:
    allowed_species_set = set(allowed_species)

    focal_matches = tree.search_nodes(name=focal_species)
    if len(focal_matches) != 1:
        raise ValueError(
            f"Expected exactly 1 tip named '{focal_species}', found {len(focal_matches)}."
        )

    current = focal_matches[0]
    branch_sets: List[Set[str]] = []
    seen = set()

    while current is not None:
        descendants = get_allowed_descendants(current, allowed_species_set)

        if focal_species not in descendants:
            raise ValueError(
                f"Focal species '{focal_species}' was not found in one of its ancestral sets."
            )

        frozen = frozenset(descendants)
        if frozen not in seen:
            seen.add(frozen)
            branch_sets.append(descendants)

        current = current.up

    return branch_sets


def write_branch_definition_file(
    focal_species: str,
    branch_sets: List[Set[str]],
    output_dir: Path,
) -> Path:
    output_path = output_dir / f"{focal_species}.branch_definitions.tsv"

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        for branch_id, species_set in enumerate(branch_sets, start=1):
            for species in sorted(species_set):
                writer.writerow([branch_id, species])

    return output_path


def main() -> None:
    args = parse_args()

    tree = read_tree(args.tree, args.tree_format)
    allowed_species = read_redip_species(
        genomes_tsv=args.genomes_tsv,
        allowed_species=args.allowed_species,
    )

    tree_tip_names = get_tree_tip_names(tree)
    missing_species = [sp for sp in allowed_species if sp not in tree_tip_names]

    if missing_species:
        raise ValueError(
            "The following rediploidisation/WGD ingroup species were not found "
            "in the species tree: "
            + ", ".join(missing_species)
        )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for focal_species in allowed_species:
        branch_sets = make_branch_sets_for_focal(
            tree=tree,
            focal_species=focal_species,
            allowed_species=allowed_species,
        )

        output_path = write_branch_definition_file(
            focal_species=focal_species,
            branch_sets=branch_sets,
            output_dir=output_dir,
        )

        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
