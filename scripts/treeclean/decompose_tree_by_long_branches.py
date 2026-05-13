#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from typing import List, Optional, Set, Tuple

from ete3 import Tree


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split a gene tree into cleaning units by repeatedly cutting the "
            "longest eligible branch. A branch is eligible if it is at least "
            "--min-branch long and both resulting sides contain at least "
            "--min-tips tips."
        )
    )
    parser.add_argument("--og", required=True, help="Original OG identifier, e.g. og_1234")
    parser.add_argument("--tree", required=True, help="Input Newick tree with branch lengths")
    parser.add_argument("--min-branch", type=float, required=True, help="Minimum branch length to cut")
    parser.add_argument("--min-tips", type=int, default=4, help="Minimum tips required on both sides of a cut")
    parser.add_argument("--tree-format", type=int, default=1, help="ETE tree format code")
    parser.add_argument("--out-membership", required=True, help="Output cleaning-unit membership TSV")
    parser.add_argument("--out-report", required=True, help="Output decomposition report TSV")
    return parser.parse_args()


def read_tree(path: str, tree_format: int) -> Tree:
    tree = Tree(path, format=tree_format)
    unnamed = [leaf for leaf in tree.iter_leaves() if not leaf.name]
    if unnamed:
        raise ValueError(f"Tree contains {len(unnamed)} unnamed terminal tips: {path}")
    return tree


def descendant_tips(node) -> Set[str]:
    return {leaf.name for leaf in node.iter_leaves()}


def find_best_split(
    tree: Tree,
    component: Set[str],
    min_branch: float,
    min_tips: int,
) -> Optional[dict]:
    best = None

    for node in tree.traverse("postorder"):
        if node.is_root():
            continue

        branch_length = float(node.dist or 0.0)
        if branch_length < min_branch:
            continue

        side_a = descendant_tips(node) & component
        if not side_a:
            continue

        side_b = component - side_a

        if len(side_a) < min_tips or len(side_b) < min_tips:
            continue

        if best is None or branch_length > best["branch_length"]:
            best = {
                "branch_length": branch_length,
                "side_a": side_a,
                "side_b": side_b,
            }

    return best


def decompose(
    tree: Tree,
    min_branch: float,
    min_tips: int,
) -> Tuple[List[Set[str]], List[dict]]:
    components = [{leaf.name for leaf in tree.iter_leaves()}]
    cuts = []

    while True:
        changed = False
        next_components = []

        for component in components:
            split = find_best_split(tree, component, min_branch, min_tips)

            if split is None:
                next_components.append(component)
                continue

            changed = True
            cuts.append(
                {
                    "cut_index": len(cuts) + 1,
                    "branch_length": split["branch_length"],
                    "n_tips_side_a": len(split["side_a"]),
                    "n_tips_side_b": len(split["side_b"]),
                }
            )
            next_components.append(split["side_a"])
            next_components.append(split["side_b"])

        components = next_components

        if not changed:
            break

    return components, cuts


def write_membership(
    path: str,
    og: str,
    components: List[Set[str]],
    was_decomposed: bool,
) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "original_og",
                "cleaning_unit_id",
                "subtree_index",
                "was_decomposed",
                "tip",
            ]
        )

        for i, tips in enumerate(components, start=1):
            if was_decomposed:
                unit = f"{og}_subtree_{i:03d}"
                subtree_index = str(i)
            else:
                unit = og
                subtree_index = "NA"

            for tip in sorted(tips):
                writer.writerow(
                    [
                        og,
                        unit,
                        subtree_index,
                        str(was_decomposed).lower(),
                        tip,
                    ]
                )


def write_report(
    path: str,
    og: str,
    n_tips_original: int,
    components: List[Set[str]],
    cuts: List[dict],
    min_branch: float,
    min_tips: int,
) -> None:
    was_decomposed = bool(cuts)

    fieldnames = [
        "original_og",
        "was_decomposed",
        "n_tips_original",
        "n_subtrees",
        "n_cuts",
        "min_branch",
        "min_tips",
        "cut_index",
        "cut_branch_length",
        "n_tips_side_a",
        "n_tips_side_b",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        if not cuts:
            writer.writerow(
                {
                    "original_og": og,
                    "was_decomposed": str(was_decomposed).lower(),
                    "n_tips_original": n_tips_original,
                    "n_subtrees": len(components),
                    "n_cuts": 0,
                    "min_branch": min_branch,
                    "min_tips": min_tips,
                    "cut_index": "NA",
                    "cut_branch_length": "NA",
                    "n_tips_side_a": "NA",
                    "n_tips_side_b": "NA",
                }
            )
            return

        for cut in cuts:
            writer.writerow(
                {
                    "original_og": og,
                    "was_decomposed": str(was_decomposed).lower(),
                    "n_tips_original": n_tips_original,
                    "n_subtrees": len(components),
                    "n_cuts": len(cuts),
                    "min_branch": min_branch,
                    "min_tips": min_tips,
                    "cut_index": cut["cut_index"],
                    "cut_branch_length": cut["branch_length"],
                    "n_tips_side_a": cut["n_tips_side_a"],
                    "n_tips_side_b": cut["n_tips_side_b"],
                }
            )


def main() -> None:
    args = parse_args()

    tree = read_tree(args.tree, args.tree_format)
    n_tips_original = len(list(tree.iter_leaves()))

    components, cuts = decompose(
        tree=tree,
        min_branch=args.min_branch,
        min_tips=args.min_tips,
    )

    was_decomposed = bool(cuts)

    write_membership(
        path=args.out_membership,
        og=args.og,
        components=components,
        was_decomposed=was_decomposed,
    )

    write_report(
        path=args.out_report,
        og=args.og,
        n_tips_original=n_tips_original,
        components=components,
        cuts=cuts,
        min_branch=args.min_branch,
        min_tips=args.min_tips,
    )


if __name__ == "__main__":
    main()
