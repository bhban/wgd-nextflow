#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect IQ-TREE treefiles from per-unit IQ-TREE directories and write "
            "a multi-tree file plus an order map for TreeShrink."
        )
    )
    parser.add_argument(
        "--iqtree-dirs",
        nargs="+",
        required=True,
        help="One or more IQ-TREE result directories",
    )
    parser.add_argument("--out-trees", required=True)
    parser.add_argument("--out-order", required=True)
    return parser.parse_args()


def find_treefile(iqtree_dir: Path) -> Path:
    candidates = sorted(iqtree_dir.glob("*_iqtree.treefile"))

    if not candidates:
        candidates = sorted(iqtree_dir.glob("*.treefile"))

    if not candidates:
        raise FileNotFoundError(f"No IQ-TREE treefile found in {iqtree_dir}")

    if len(candidates) > 1:
        preferred = [p for p in candidates if p.name == f"{iqtree_dir.name}_iqtree.treefile"]
        if len(preferred) == 1:
            return preferred[0]

    return candidates[0]


def read_newick(path: Path) -> str:
    text = path.read_text().strip()
    if not text:
        raise ValueError(f"Empty treefile: {path}")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) != 1:
        return "".join(lines)

    return lines[0]


def main() -> None:
    args = parse_args()

    records: List[Tuple[int, str, Path, str]] = []

    for i, d in enumerate(sorted({str(Path(x)) for x in args.iqtree_dirs}), start=1):
        iqtree_dir = Path(d)
        treefile = find_treefile(iqtree_dir)
        unit = iqtree_dir.name
        newick = read_newick(treefile)
        records.append((i, unit, treefile, newick))

    with open(args.out_trees, "w") as handle:
        for _, _, _, newick in records:
            handle.write(newick.rstrip() + "\n")

    with open(args.out_order, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tree_index", "cleaning_unit_id", "treefile"])
        for tree_index, unit, treefile, _ in records:
            writer.writerow([tree_index, unit, treefile])


if __name__ == "__main__":
    main()
