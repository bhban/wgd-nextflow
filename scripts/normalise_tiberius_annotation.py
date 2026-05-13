#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


def format_attrs(attrs: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in attrs.items())


def transcript_to_gene(transcript_id: str) -> str:
    match = re.match(r"(.+)\.t\d+$", transcript_id)
    if match:
        return match.group(1)
    return transcript_id


def make_region_id(region_id: str, feature_id: str) -> str:
    """
    Make an ID unique within the genome using the GFF column-1 region/chromosome ID.

    Example:
        region_id = CDWWSV010000001.1
        feature_id = g1.t1
        output = CDWWSV010000001.1_g1.t1
    """
    prefix = f"{region_id}_"
    if feature_id.startswith(prefix):
        return feature_id
    return f"{prefix}{feature_id}"


def parse_tiberius_fasta_header(header: str) -> tuple[str, str]:
    """
    Tiberius FASTA headers look like:

        >g1|g1.t1|CDWWSV010000001.1:4371-7508(-)

    Return:
        region_id = CDWWSV010000001.1
        transcript_id = g1.t1
    """
    first = header.strip().split()[0]
    if first.startswith(">"):
        first = first[1:]

    parts = first.split("|")

    if len(parts) >= 3:
        transcript_id = parts[1]
        region_field = parts[2]
        region_id = region_field.split(":", 1)[0]
        return region_id, transcript_id

    if len(parts) >= 2:
        return "unknown_region", parts[1]

    return "unknown_region", first


def normalise_fasta(in_fasta: Path, out_fasta: Path) -> None:
    with in_fasta.open() as inp, out_fasta.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                region_id, transcript_id = parse_tiberius_fasta_header(line)
                out.write(f">{make_region_id(region_id, transcript_id)}\n")
            else:
                out.write(line.upper())


def normalise_gff3(in_gff: Path, out_gff: Path) -> None:
    transcript_map: dict[tuple[str, str], str] = {}

    with in_gff.open() as inp, out_gff.open("w") as out:
        out.write("##gff-version 3\n")

        for raw_line in inp:
            line = raw_line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attr_text = fields
            attrs = parse_attrs(attr_text)

            if feature == "mRNA":
                old_tx = attrs.get("ID")

                if old_tx is None:
                    raise ValueError(f"mRNA feature has no ID: {line}")

                old_gene = transcript_to_gene(old_tx)

                new_gene = make_region_id(seqid, old_gene)
                new_tx = make_region_id(seqid, old_tx)

                transcript_map[(seqid, old_tx)] = new_tx

                gene_fields = [
                    seqid,
                    source,
                    "gene",
                    start,
                    end,
                    score,
                    strand,
                    ".",
                    format_attrs({"ID": new_gene, "Name": new_gene}),
                ]
                out.write("\t".join(gene_fields) + "\n")

                fields[8] = format_attrs(
                    {
                        "ID": new_tx,
                        "Parent": new_gene,
                        "Name": new_tx,
                    }
                )
                out.write("\t".join(fields) + "\n")

            elif feature in {"CDS", "exon"}:
                old_parent = attrs.get("Parent")

                if old_parent is None:
                    continue

                new_parent = transcript_map.get(
                    (seqid, old_parent),
                    make_region_id(seqid, old_parent),
                )

                old_id = attrs.get("ID")

                if old_id:
                    new_id = make_region_id(seqid, old_id)
                else:
                    new_id = f"{new_parent}.{feature}"

                fields[8] = format_attrs(
                    {
                        "ID": new_id,
                        "Parent": new_parent,
                    }
                )
                out.write("\t".join(fields) + "\n")

            else:
                if "ID" in attrs:
                    attrs["ID"] = make_region_id(seqid, attrs["ID"])

                if "Parent" in attrs:
                    attrs["Parent"] = transcript_map.get(
                        (seqid, attrs["Parent"]),
                        make_region_id(seqid, attrs["Parent"]),
                    )

                fields[8] = format_attrs(attrs)
                out.write("\t".join(fields) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalise Tiberius GFF3/CDS/PEP outputs using GFF column-1 region IDs."
    )
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--cds", required=True, type=Path)
    parser.add_argument("--pep", required=True, type=Path)
    parser.add_argument("--out-gff", required=True, type=Path)
    parser.add_argument("--out-cds", required=True, type=Path)
    parser.add_argument("--out-pep", required=True, type=Path)
    args = parser.parse_args()

    args.out_gff.parent.mkdir(parents=True, exist_ok=True)
    args.out_cds.parent.mkdir(parents=True, exist_ok=True)
    args.out_pep.parent.mkdir(parents=True, exist_ok=True)

    normalise_gff3(args.gff, args.out_gff)
    normalise_fasta(args.cds, args.out_cds)
    normalise_fasta(args.pep, args.out_pep)


if __name__ == "__main__":
    main()
