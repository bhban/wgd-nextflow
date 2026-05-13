#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path
from typing import Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select FASTAs whose cleaning_unit_id is present in a passing membership file."
    )
    parser.add_argument("--fastas", nargs="+", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-manifest", required=True)
    return parser.parse_args()


def read_passing_units(path: str) -> Set[str]:
    units = set()

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            return units

        if "cleaning_unit_id" not in reader.fieldnames:
            raise ValueError("Membership file must contain cleaning_unit_id")

        for row in reader:
            unit = row["cleaning_unit_id"].strip()
            if unit:
                units.add(unit)

    return units


def unit_from_fasta(path: Path) -> str:
    name = path.name

    if name.endswith(".pruned.fasta"):
        return name.removesuffix(".pruned.fasta")

    if name.endswith(".fasta"):
        return name.removesuffix(".fasta")

    return path.stem


def main() -> None:
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    passing_units = read_passing_units(args.membership)

    with open(args.out_manifest, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cleaning_unit_id", "source_fasta", "selected_fasta", "selected"])

        for fasta in sorted(args.fastas):
            fasta_path = Path(fasta)
            unit = unit_from_fasta(fasta_path)

            if unit not in passing_units:
                writer.writerow([unit, fasta_path, "NA", "false"])
                continue

            out_fasta = out_dir / f"{unit}.fasta"
            shutil.copyfile(fasta_path, out_fasta)
            writer.writerow([unit, fasta_path, out_fasta, "true"])


if __name__ == "__main__":
    main()
