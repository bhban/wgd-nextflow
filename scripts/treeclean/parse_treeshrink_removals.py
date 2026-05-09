#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse TreeShrink removal/outlier outputs into a simple TSV with "
            "cleaning_unit_id and tip. If no recognizable removal file is found, "
            "an empty removals TSV is written."
        )
    )
    parser.add_argument("--treeshrink-dir", required=True)
    parser.add_argument("--prefix", default="treeshrink")
    parser.add_argument("--alpha", default="0.005")
    parser.add_argument("--tree-order", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def read_tree_order(path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    index_to_unit = {}
    unit_to_treefile = {}

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"tree_index", "cleaning_unit_id", "treefile"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Tree order file missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            index_to_unit[str(row["tree_index"])] = row["cleaning_unit_id"]
            unit_to_treefile[row["cleaning_unit_id"]] = row["treefile"]

    return index_to_unit, unit_to_treefile


def tokenize(line: str) -> List[str]:
    cleaned = (
        line.strip()
        .replace(",", " ")
        .replace(";", " ")
        .replace("[", " ")
        .replace("]", " ")
        .replace("(", " ")
        .replace(")", " ")
        .replace("{", " ")
        .replace("}", " ")
        .replace(":", " ")
    )
    return [tok.strip() for tok in re.split(r"\s+", cleaned) if tok.strip()]


def likely_removal_file(path: Path) -> bool:
    name = path.name.lower()
    if path.is_dir():
        return False
    if path.stat().st_size == 0:
        return False
    if any(x in name for x in ["remov", "outlier", "shrink", "output"]):
        return True
    return False


def candidate_files(root: Path) -> List[Path]:
    files = []
    for path in root.rglob("*"):
        if likely_removal_file(path):
            files.append(path)
    return sorted(files)


def parse_table_like_file(
    path: Path,
    index_to_unit: Dict[str, str],
    known_units: Set[str],
) -> List[Tuple[str, str, str]]:
    removals = []

    try:
        text = path.read_text(errors="replace")
    except UnicodeDecodeError:
        return removals

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        lower = line.lower()
        if lower.startswith(("tree", "gene", "family", "unit")) and "remove" in lower:
            continue

        tokens = tokenize(line)
        if len(tokens) < 2:
            continue

        first = tokens[0]
        unit = None

        if first in known_units:
            unit = first
            tips = tokens[1:]
        elif first in index_to_unit:
            unit = index_to_unit[first]
            tips = tokens[1:]
        else:
            tips = tokens

            for tok in tokens:
                if tok in known_units:
                    unit = tok
                    tips = [x for x in tokens if x != tok]
                    break

        if unit is None:
            continue

        for tip in tips:
            if tip in {"NA", "None", "none", "null", "-", "remove", "removed", "tips", "taxa"}:
                continue

            if tip in known_units or tip in index_to_unit:
                continue

            if re.fullmatch(r"[0-9.]+", tip):
                continue

            removals.append((unit, tip, str(path)))

    return removals


def main() -> None:
    args = parse_args()

    treeshrink_dir = Path(args.treeshrink_dir)
    index_to_unit, unit_to_treefile = read_tree_order(args.tree_order)
    known_units = set(unit_to_treefile)

    all_removals = []

    for path in candidate_files(treeshrink_dir):
        all_removals.extend(parse_table_like_file(path, index_to_unit, known_units))

    unique = sorted(set(all_removals))

    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cleaning_unit_id", "tip", "source_file"])
        for unit, tip, source_file in unique:
            writer.writerow([unit, tip, source_file])


if __name__ == "__main__":
    main()
