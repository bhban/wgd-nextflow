#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from typing import Dict, List, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter tree-cleaning units by species count and outgroup presence. "
            "The input membership TSV must contain original_og, cleaning_unit_id, "
            "subtree_index, and tip."
        )
    )
    parser.add_argument("--membership", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--min-species", type=int, default=4)
    parser.add_argument("--require-outgroup", action="store_true")
    parser.add_argument("--tip-separator", default="|")
    parser.add_argument(
        "--tip-label-format",
        default="species_chr_gene",
        choices=["species_chr_gene", "species_gene", "species_only", "full_label"],
    )
    parser.add_argument("--out-membership", required=True)
    parser.add_argument("--out-report", required=True)
    return parser.parse_args()


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "outgroup"}


def parse_species(tip: str, separator: str, label_format: str) -> str:
    if label_format in {"species_chr_gene", "species_gene"}:
        return tip.split(separator)[0]
    if label_format == "species_only":
        return tip
    if label_format == "full_label":
        return tip
    raise ValueError(f"Unsupported tip label format: {label_format}")


def read_outgroups(genomes_tsv: str) -> Set[str]:
    outgroups = set()

    with open(genomes_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Could not read header from genomes TSV: {genomes_tsv}")

        if "genome_id" not in reader.fieldnames:
            raise ValueError("genomes.tsv must contain a genome_id column")

        has_outgroup_col = "outgroup" in reader.fieldnames

        for row in reader:
            genome = row["genome_id"].strip()
            if not genome:
                continue

            if has_outgroup_col and truthy(row.get("outgroup", "0")):
                outgroups.add(genome)

    return outgroups


def main() -> None:
    args = parse_args()
    outgroups = read_outgroups(args.genomes_tsv)

    rows_by_unit: Dict[str, List[dict]] = defaultdict(list)
    tips_by_unit: Dict[str, Set[str]] = defaultdict(set)
    original_og_by_unit: Dict[str, str] = {}
    subtree_index_by_unit: Dict[str, str] = {}

    with open(args.membership, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"original_og", "cleaning_unit_id", "subtree_index", "tip"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Membership file missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            unit = row["cleaning_unit_id"]
            tip = row["tip"]
            rows_by_unit[unit].append(row)
            tips_by_unit[unit].add(tip)
            original_og_by_unit[unit] = row["original_og"]
            subtree_index_by_unit[unit] = row["subtree_index"]

    passing_units = set()
    report_rows = []

    for unit in sorted(tips_by_unit):
        tips = tips_by_unit[unit]
        species = {
            parse_species(tip, args.tip_separator, args.tip_label_format)
            for tip in tips
        }
        present_outgroups = sorted(species & outgroups)

        fail_reasons = []

        if len(species) < args.min_species:
            fail_reasons.append("too_few_species")

        if args.require_outgroup and not present_outgroups:
            fail_reasons.append("no_outgroup")

        passes = not fail_reasons

        if passes:
            passing_units.add(unit)

        report_rows.append(
            {
                "original_og": original_og_by_unit[unit],
                "cleaning_unit_id": unit,
                "subtree_index": subtree_index_by_unit[unit],
                "n_tips": len(tips),
                "n_species": len(species),
                "species": ",".join(sorted(species)) if species else "NA",
                "has_outgroup": str(bool(present_outgroups)).lower(),
                "outgroup_species": ",".join(present_outgroups) if present_outgroups else "NA",
                "passes_filter": str(passes).lower(),
                "fail_reason": ";".join(fail_reasons) if fail_reasons else "NA",
            }
        )

    with open(args.out_membership, "w", newline="") as handle:
        fieldnames = ["original_og", "cleaning_unit_id", "subtree_index", "tip"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for unit in sorted(passing_units):
            for row in rows_by_unit[unit]:
                writer.writerow({key: row[key] for key in fieldnames})

    with open(args.out_report, "w", newline="") as handle:
        fieldnames = [
            "original_og",
            "cleaning_unit_id",
            "subtree_index",
            "n_tips",
            "n_species",
            "species",
            "has_outgroup",
            "outgroup_species",
            "passes_filter",
            "fail_reason",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)


if __name__ == "__main__":
    main()
