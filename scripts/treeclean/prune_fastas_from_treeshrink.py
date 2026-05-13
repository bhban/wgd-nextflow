#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove TreeShrink-flagged tips from cleaning-unit raw FASTAs."
    )
    parser.add_argument("--fastas", nargs="+", required=True)
    parser.add_argument("--removals", required=True)
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--tip-separator", default="|")
    parser.add_argument(
        "--tip-label-format",
        default="species_chr_gene",
        choices=["species_chr_gene", "species_gene", "species_only", "full_label"],
    )
    parser.add_argument("--protect-outgroups", action="store_true")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-membership", required=True)
    parser.add_argument("--out-report", required=True)
    return parser.parse_args()


def truthy(value: str) -> bool:
    text = str(value).strip().lower()

    if text in {"", "0", "false", "f", "no", "n", "none", "na", "nan"}:
        return False

    if text in {"true", "t", "yes", "y", "outgroup"}:
        return True

    try:
        return float(text) > 0
    except ValueError:
        return False


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


def read_fasta(path: str) -> Dict[str, str]:
    seqs = {}
    current = None
    chunks: List[str] = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if current is not None:
        seqs[current] = "".join(chunks)

    return seqs


def write_fasta(path: Path, seqs: Dict[str, str]) -> None:
    with open(path, "w") as handle:
        for name in sorted(seqs):
            handle.write(f">{name}\n")
            seq = seqs[name]
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def read_removals(path: str) -> Dict[str, Set[str]]:
    removals = defaultdict(set)

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            return removals

        required = {"cleaning_unit_id", "tip"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Removal file missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            unit = row["cleaning_unit_id"].strip()
            tip = row["tip"].strip()
            if unit and tip:
                removals[unit].add(tip)

    return removals


def parse_original_og(unit: str) -> str:
    marker = "_subtree_"
    if marker in unit:
        return unit.split(marker)[0]
    return unit


def parse_subtree_index(unit: str) -> str:
    marker = "_subtree_"
    if marker in unit:
        return unit.split(marker, 1)[1]
    return "NA"


def was_decomposed(unit: str) -> str:
    return str("_subtree_" in unit).lower()


def main() -> None:
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    removals = read_removals(args.removals)
    outgroups = read_outgroups(args.genomes_tsv)

    report_rows = []

    with open(args.out_membership, "w", newline="") as membership_handle:
        membership_writer = csv.writer(membership_handle, delimiter="\t")
        membership_writer.writerow(
            [
                "original_og",
                "cleaning_unit_id",
                "subtree_index",
                "was_decomposed",
                "was_pruned",
                "tip",
            ]
        )

        for fasta in sorted(args.fastas):
            fasta_path = Path(fasta)
            unit = fasta_path.stem
            seqs = read_fasta(fasta)

            tips_before = set(seqs)
            requested_removals = removals.get(unit, set())

            protected_requested = set()
            if args.protect_outgroups:
                for tip in requested_removals:
                    species = parse_species(tip, args.tip_separator, args.tip_label_format)
                    if species in outgroups:
                        protected_requested.add(tip)

            effective_requested_removals = requested_removals - protected_requested

            actual_removals = tips_before & effective_requested_removals
            missing_requested = effective_requested_removals - tips_before
            kept = {tip: seq for tip, seq in seqs.items() if tip not in actual_removals}

            was_pruned = bool(actual_removals)

            out_fasta = out_dir / f"{unit}.pruned.fasta"
            write_fasta(out_fasta, kept)

            original_og = parse_original_og(unit)
            subtree_index = parse_subtree_index(unit)
            decomposed = was_decomposed(unit)

            for tip in sorted(kept):
                membership_writer.writerow(
                    [
                        original_og,
                        unit,
                        subtree_index,
                        decomposed,
                        str(was_pruned).lower(),
                        tip,
                    ]
                )

            report_rows.append(
                {
                    "original_og": original_og,
                    "cleaning_unit_id": unit,
                    "subtree_index": subtree_index,
                    "was_decomposed": decomposed,
                    "was_pruned": str(was_pruned).lower(),
                    "fasta": str(fasta_path),
                    "pruned_fasta": str(out_fasta),
                    "n_tips_before": len(tips_before),
                    "n_requested_removals": len(requested_removals),
                    "n_protected_requested_removals": len(protected_requested),
                    "n_effective_requested_removals": len(effective_requested_removals),
                    "n_removed": len(actual_removals),
                    "n_missing_requested_removals": len(missing_requested),
                    "n_tips_after": len(kept),
                    "removed_tips": ",".join(sorted(actual_removals)) if actual_removals else "NA",
                    "protected_requested_removals": ",".join(sorted(protected_requested)) if protected_requested else "NA",
                    "missing_requested_removals": ",".join(sorted(missing_requested)) if missing_requested else "NA",
                }
            )

    with open(args.out_report, "w", newline="") as handle:
        fieldnames = [
            "original_og",
            "cleaning_unit_id",
            "subtree_index",
            "was_decomposed",
            "was_pruned",
            "fasta",
            "pruned_fasta",
            "n_tips_before",
            "n_requested_removals",
            "n_protected_requested_removals",
            "n_effective_requested_removals",
            "n_removed",
            "n_missing_requested_removals",
            "n_tips_after",
            "removed_tips",
            "protected_requested_removals",
            "missing_requested_removals",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)


if __name__ == "__main__":
    main()
