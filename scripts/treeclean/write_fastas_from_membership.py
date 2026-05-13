#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write one raw FASTA per passing cleaning unit from a membership TSV."
    )
    parser.add_argument("--og", required=True)
    parser.add_argument("--raw-fasta", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-manifest", required=True)
    return parser.parse_args()


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
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


def write_fasta(path: Path, records: Dict[str, str]) -> None:
    with open(path, "w") as handle:
        for name in sorted(records):
            handle.write(f">{name}\n")
            seq = records[name]
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main() -> None:
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = read_fasta(args.raw_fasta)

    tips_by_unit = defaultdict(list)
    original_og_by_unit = {}
    subtree_index_by_unit = {}
    was_decomposed_by_unit = {}

    with open(args.membership, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"original_og", "cleaning_unit_id", "subtree_index", "tip"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Membership file missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            unit = row["cleaning_unit_id"]
            tips_by_unit[unit].append(row["tip"])
            original_og_by_unit[unit] = row["original_og"]
            subtree_index_by_unit[unit] = row["subtree_index"]
            was_decomposed_by_unit[unit] = row.get("was_decomposed", "NA")

    with open(args.out_manifest, "w", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t")
        writer.writerow(
            [
                "original_og",
                "cleaning_unit_id",
                "subtree_index",
                "was_decomposed",
                "fasta",
                "n_tips",
                "status",
                "message",
            ]
        )

        for unit in sorted(tips_by_unit):
            tips = sorted(set(tips_by_unit[unit]))
            missing_tips = [tip for tip in tips if tip not in seqs]

            if missing_tips:
                writer.writerow(
                    [
                        original_og_by_unit[unit],
                        unit,
                        subtree_index_by_unit[unit],
                        was_decomposed_by_unit[unit],
                        "NA",
                        len(tips),
                        "FAIL",
                        f"missing_tips:{','.join(missing_tips)}",
                    ]
                )
                raise KeyError(
                    f"{len(missing_tips)} tips from {unit} were not found in {args.raw_fasta}: "
                    + ", ".join(missing_tips[:10])
                )

            out_fasta = out_dir / f"{unit}.fasta"
            records = {tip: seqs[tip] for tip in tips}
            write_fasta(out_fasta, records)

            writer.writerow(
                [
                    original_og_by_unit[unit],
                    unit,
                    subtree_index_by_unit[unit],
                    was_decomposed_by_unit[unit],
                    out_fasta,
                    len(tips),
                    "OK",
                    "NA",
                ]
            )


if __name__ == "__main__":
    main()
