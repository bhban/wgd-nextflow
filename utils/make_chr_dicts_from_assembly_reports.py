#!/usr/bin/env python3

import argparse
import csv
import gzip
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def clean_chrom_name(value, prefix="chr"):
    value = value.strip()

    if value.lower() in {"na", "none", ""}:
        return None

    if value.lower().startswith("chr"):
        return "chr" + value[3:]

    return f"{prefix}{value}"


def is_chromosome_row(sequence_role, assigned_type, assembly_unit):
    text = " ".join([
        sequence_role or "",
        assigned_type or "",
        assembly_unit or "",
    ]).lower()

    return "chromosome" in text or sequence_role == "assembled-molecule"


def make_chr_dict(report_path, output_path, chrom_prefix="chr"):
    rows_written = 0

    with open_maybe_gzip(report_path) as handle, open(output_path, "w") as out:
        out.write("#chrom\tchromStart\tchromEnd\tname\n")

        for line in handle:
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 10:
                continue

            sequence_role = fields[1]
            assigned_molecule = fields[2]
            assigned_type = fields[3]
            genbank_accn = fields[4]
            refseq_accn = fields[6]
            assembly_unit = fields[7]
            sequence_length = fields[8]

            if genbank_accn.lower() == "na" and refseq_accn.lower() == "na":
                continue

            try:
                chrom_end = int(sequence_length)
            except ValueError:
                continue

            # Use the RefSeq accession when available, because RefSeq GFFs use
            # RefSeq sequence IDs such as NC_064543.1 rather than GenBank IDs
            # such as CM043645.1. Fall back to GenBank when no RefSeq accession
            # is present.
            old_seqid = (
                refseq_accn
                if refseq_accn.lower() != "na"
                else genbank_accn
            )

            if is_chromosome_row(sequence_role, assigned_type, assembly_unit):
                new_seqid = clean_chrom_name(assigned_molecule, prefix=chrom_prefix)
                if new_seqid is None:
                    new_seqid = old_seqid
            else:
                new_seqid = old_seqid

            out.write(f"{new_seqid}\t0\t{chrom_end}\t{old_seqid}\n")
            rows_written += 1

    return rows_written


def download_file(url, output_path):
    with urllib.request.urlopen(url) as response, open(output_path, "wb") as out:
        out.write(response.read())


def read_manifest(manifest_path):
    with open(manifest_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"genome_id", "assembly_report_url"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"Manifest is missing required column(s): {', '.join(sorted(missing))}"
            )

        for row in reader:
            genome_id = row["genome_id"].strip()
            url = row["assembly_report_url"].strip()

            if not genome_id or not url:
                continue

            yield genome_id, url


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download NCBI assembly reports from a TSV manifest and write "
            "Style B chr_dict files: #chrom chromStart chromEnd name."
        )
    )

    parser.add_argument(
        "manifest",
        help="TSV with columns: genome_id, assembly_report_url",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="chr_dict",
        help="Output directory for chr_dict files. Default: chr_dict",
    )

    parser.add_argument(
        "--reports-dir",
        default=None,
        help=(
            "Optional directory to save downloaded assembly reports. "
            "If omitted, reports are downloaded to a temporary directory."
        ),
    )

    parser.add_argument(
        "--chrom-prefix",
        default="chr",
        help="Prefix for chromosome names when rows are chromosome-level. Default: chr",
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.reports_dir:
        reports_dir = Path(args.reports_dir)
        reports_dir.mkdir(parents=True, exist_ok=True)
        temp_context = None
    else:
        temp_context = tempfile.TemporaryDirectory()
        reports_dir = Path(temp_context.name)

    n_done = 0
    n_fail = 0

    try:
        for genome_id, url in read_manifest(args.manifest):
            report_path = reports_dir / f"{genome_id}_assembly_report.txt"
            output_path = outdir / f"{genome_id}.tsv"

            print(f"Downloading assembly report for {genome_id}", file=sys.stderr)

            try:
                download_file(url, report_path)
                rows_written = make_chr_dict(
                    report_path=report_path,
                    output_path=output_path,
                    chrom_prefix=args.chrom_prefix,
                )

                print(
                    f"  wrote {output_path} with {rows_written} rows",
                    file=sys.stderr,
                )
                n_done += 1

            except Exception as error:
                print(
                    f"  failed for {genome_id}: {error}",
                    file=sys.stderr,
                )
                n_fail += 1

    finally:
        if temp_context is not None:
            temp_context.cleanup()

    print(
        f"Finished: {n_done} chr_dict files written; {n_fail} failed",
        file=sys.stderr,
    )

    if n_fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
