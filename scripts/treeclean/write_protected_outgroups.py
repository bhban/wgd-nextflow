#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write comma/newline compatible list of outgroup species from genomes.tsv."
    )
    parser.add_argument("--genomes-tsv", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "outgroup"}


def main() -> None:
    args = parse_args()

    outgroups = []

    with open(args.genomes_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise ValueError(f"Could not read header from genomes TSV: {args.genomes_tsv}")

        if "genome_id" not in reader.fieldnames:
            raise ValueError("genomes.tsv must contain a genome_id column")

        if "outgroup" not in reader.fieldnames:
            outgroups = []
        else:
            for row in reader:
                genome = row["genome_id"].strip()
                if genome and truthy(row.get("outgroup", "0")):
                    outgroups.append(genome)

    with open(args.out, "w") as handle:
        for genome in sorted(set(outgroups)):
            handle.write(genome + "\n")


if __name__ == "__main__":
    main()
