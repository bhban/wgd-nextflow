#!/usr/bin/env python3

"""
Adapted from Classify_oncorhynchus_mykiss_rediploidization_histories_2026_03_03.py
(Written by Drew Larson)

Read the classification TSV produced by classify_redip_events.py, classify each
row into the smallest matching rediploidization branch, look up genomic positions
for the two target genes, and write a Circos links file.

Supports two tip-label styles:
    species@gene_id
    species|chr|gene_id

Branch definitions are supplied in a TSV file with two columns:
    branch_id    species

Example branch-definitions TSV:
    1   Oncorhynchus_mykiss
    2   Oncorhynchus_mykiss
    2   Oncorhynchus_kisutch
    3   Oncorhynchus_mykiss
    3   Oncorhynchus_kisutch
    3   Salvelinus_alpinus
    3   Salvelinus_leucomaenis

Optional colors TSV format:
    branch_id    color

Example colors TSV:
    1   color=245,235,39
    2   color=248,149,64
    3   color=204,71,120

Position tables can use arbitrary column names, provided you specify them with:
    --position-key-column
    --position-chr-column
    --position-start-column
    --position-end-column

For example, a table with columns:
    ofID    genome    og    flag    id    chr    start    end    ord

can be used with:
    --position-key-column id
    --position-chr-column chr
    --position-start-column start
    --position-end-column end
    --position-species-column genome

Key is usually the gene_id, but can also be the full tip label if you
set --position-key-type full_label.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


# ----------------------------------------------------------------------
# Argument parsing and logging
# ----------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert rediploidization classifications into a Circos links file."
        )
    )

    parser.add_argument(
        "--classification-tsv",
        required=True,
        help="Classification TSV produced by classify_redip_events.py."
    )
    parser.add_argument(
        "--positions",
        required=True,
        help="Gene-position TSV."
    )
    parser.add_argument(
        "--branch-definitions",
        required=True,
        help="TSV defining branch_id to species memberships."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV for Circos links."
    )
    parser.add_argument(
        "--colors",
        default=None,
        help="Optional TSV mapping branch_id to color."
    )
    parser.add_argument(
        "--tip-separator",
        default="@",
        help="Separator used in target tip labels (default: @)."
    )
    parser.add_argument(
        "--label-format",
        choices=["species_gene", "species_chr_gene"],
        default="species_gene",
        help="Tip-label format (default: species_gene)."
    )
    parser.add_argument(
        "--position-key-type",
        choices=["gene", "full_label"],
        default="gene",
        help=(
            "How to match target tips to the positions table: "
            "'gene' uses just the parsed gene_id, "
            "'full_label' uses the whole tip label (default: gene)."
        )
    )
    parser.add_argument(
        "--position-has-header",
        action="store_true",
        help="Indicate that the positions table has a header row."
    )
    parser.add_argument(
        "--position-key-column",
        default="id",
        help="Column name in the positions file to use as the lookup key (default: id)."
    )
    parser.add_argument(
        "--position-chr-column",
        default="chr",
        help="Column name in the positions file for chromosome (default: chr)."
    )
    parser.add_argument(
        "--position-start-column",
        default="start",
        help="Column name in the positions file for start coordinate (default: start)."
    )
    parser.add_argument(
        "--position-end-column",
        default="end",
        help="Column name in the positions file for end coordinate (default: end)."
    )
    parser.add_argument(
        "--position-species-column",
        default=None,
        help=(
            "Optional column name in the positions file for species/genome. "
            "Useful for validating that the position row matches the tip label species."
        )
    )
    parser.add_argument(
        "--classification-no-header",
        action="store_true",
        help="Indicate that the classification TSV does not have a header row."
    )
    parser.add_argument(
        "--write-header",
        action="store_true",
        help="Write a header row to the output file."
    )
    parser.add_argument(
        "--include-metadata",
        action="store_true",
        help="Append metadata columns after the Circos link columns."
    )
    parser.add_argument(
        "--on-exists",
        choices=["overwrite", "append", "error"],
        default="overwrite",
        help="Behavior if output file already exists (default: overwrite)."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit immediately on malformed input instead of warning and skipping."
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level (default: INFO)."
    )

    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s: %(message)s"
    )


def handle_error(message: str, strict: bool) -> None:
    if strict:
        raise ValueError(message)
    logging.warning(message)


# ----------------------------------------------------------------------
# Parsing helpers
# ----------------------------------------------------------------------

def parse_tip_label(label: str, separator: str, label_format: str) -> Dict[str, str | None]:
    """
    Parse a tip label into its components.

    Supported formats:
        species_gene      -> species@gene_id
        species_chr_gene  -> species|chr|gene_id
    """
    parts = label.split(separator)

    if label_format == "species_gene":
        if len(parts) != 2:
            raise ValueError(
                f"Tip label '{label}' does not match format 'species_gene' "
                f"with separator '{separator}'."
            )
        species, gene_id = parts
        if not species or not gene_id:
            raise ValueError(f"Tip label '{label}' contains an empty field.")
        return {"species": species, "chr": None, "gene": gene_id}

    if label_format == "species_chr_gene":
        if len(parts) != 3:
            raise ValueError(
                f"Tip label '{label}' does not match format 'species_chr_gene' "
                f"with separator '{separator}'."
            )
        species, chr_id, gene_id = parts
        if not species or not chr_id or not gene_id:
            raise ValueError(f"Tip label '{label}' contains an empty field.")
        return {"species": species, "chr": chr_id, "gene": gene_id}

    raise ValueError(f"Unknown label format: {label_format}")


def branch_sort_key(branch_id: str) -> Tuple[int, str]:
    """
    Sort branch IDs numerically when possible, otherwise lexicographically.
    """
    try:
        return (0, f"{int(branch_id):010d}")
    except ValueError:
        return (1, branch_id)


def extract_position_key(
    tip_label: str,
    tip_separator: str,
    label_format: str,
    position_key_type: str
) -> str:
    if position_key_type == "full_label":
        return tip_label

    parsed = parse_tip_label(tip_label, tip_separator, label_format)
    gene_id = parsed["gene"]
    if gene_id is None:
        raise ValueError(f"Could not extract gene ID from tip label '{tip_label}'.")
    return gene_id


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------

def load_branch_definitions(path: str) -> Dict[str, set[str]]:
    """
    Load branch definitions from a TSV with:
        branch_id    species
    """
    branch_to_species: Dict[str, set[str]] = {}

    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(
                    f"Branch-definitions file '{path}' has fewer than 2 columns "
                    f"on line {line_number}."
                )

            branch_id = row[0].strip()
            species = row[1].strip()

            if not branch_id or not species:
                raise ValueError(
                    f"Branch-definitions file '{path}' has an empty field on line {line_number}."
                )

            branch_to_species.setdefault(branch_id, set()).add(species)

    if not branch_to_species:
        raise ValueError(f"No branch definitions were loaded from: {path}")

    return branch_to_species


def load_colors(path: str | None, branch_ids: Sequence[str]) -> Dict[str, str]:
    """
    Load colors from a TSV with:
        branch_id    color

    If no colors file is supplied, assign colors from a colourblind-friendly
    ordered palette. When only a subset of branches is present, spread those
    branches across the full palette so they remain distinguishable.
    """
    if path is None:
        full_palette = [
            "color=68,1,84",      # deep purple
            "color=72,40,120",    # purple
            "color=62,73,137",    # indigo
            "color=49,104,142",   # blue
            "color=38,130,142",   # blue-teal
            "color=31,158,137",   # teal
            "color=53,183,121",   # green
            "color=109,205,89",   # yellow-green
            "color=180,222,44",   # lime
            "color=253,231,37",   # yellow
        ]

        sorted_branch_ids = sorted(branch_ids, key=branch_sort_key)
        n_branches = len(sorted_branch_ids)
        n_colors = len(full_palette)

        if n_branches > n_colors:
            raise ValueError(
                f"No --colors file was supplied, but there are {n_branches} branch IDs "
                f"and only {n_colors} default colors."
            )

        if n_branches == 1:
            palette_indices = [n_colors // 2]
        else:
            palette_indices = [
                round(i * (n_colors - 1) / (n_branches - 1))
                for i in range(n_branches)
            ]

        return {
            branch_id: full_palette[idx]
            for branch_id, idx in zip(sorted_branch_ids, palette_indices)
        }

    color_dict: Dict[str, str] = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(
                    f"Colors file '{path}' has fewer than 2 columns on line {line_number}."
                )

            branch_id = row[0].strip()
            color = row[1].strip()

            if not branch_id or not color:
                raise ValueError(
                    f"Colors file '{path}' has an empty field on line {line_number}."
                )

            color_dict[branch_id] = color

    missing = [branch_id for branch_id in branch_ids if branch_id not in color_dict]
    if missing:
        raise ValueError(
            "Colors file does not define colors for these branch IDs: "
            f"{', '.join(sorted(missing, key=branch_sort_key))}"
        )

    return {branch_id: color_dict[branch_id] for branch_id in branch_ids}


def load_positions(
    path: str,
    has_header: bool,
    key_column: str,
    chr_column: str,
    start_column: str,
    end_column: str,
    species_column: str | None = None
) -> Dict[str, Dict[str, str]]:
    """
    Load positions from a TSV.

    Supports arbitrary column names, for example:
        ofID    genome    og    flag    id    chr    start    end    ord

    Example mapping:
        key_column="id"
        chr_column="chr"
        start_column="start"
        end_column="end"
        species_column="genome"

    Returns:
        {
            key: {
                "chr": ...,
                "start": ...,
                "end": ...,
                "species": ... or ""
            }
        }

    Duplicate keys are allowed only if all position information matches exactly.
    If a duplicate key has conflicting position information, raise an error.
    """
    pos_dict: Dict[str, Dict[str, str]] = {}

    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t") if has_header else None

        if has_header:
            if reader is None or reader.fieldnames is None:
                raise ValueError(f"Positions file '{path}' appears to have no header.")

            required = [key_column, chr_column, start_column, end_column]
            if species_column is not None:
                required.append(species_column)

            missing = [col for col in required if col not in reader.fieldnames]
            if missing:
                raise ValueError(
                    f"Positions file '{path}' is missing required columns: {', '.join(missing)}"
                )

            for line_number, row in enumerate(reader, start=2):
                key = (row.get(key_column) or "").strip()
                chrom = (row.get(chr_column) or "").strip()
                start = (row.get(start_column) or "").strip()
                end = (row.get(end_column) or "").strip()
                species = (row.get(species_column) or "").strip() if species_column else ""

                if not key or not chrom or not start or not end:
                    raise ValueError(
                        f"Positions file '{path}' has an empty required field on line {line_number}."
                    )

                new_value = {
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "species": species,
                }

                if key in pos_dict:
                    if pos_dict[key] != new_value:
                        raise ValueError(
                            f"Conflicting position information for key '{key}' in positions "
                            f"file '{path}' on line {line_number}. First value: {pos_dict[key]}; "
                            f"new value: {new_value}"
                        )
                else:
                    pos_dict[key] = new_value

        else:
            raw_reader = csv.reader(handle, delimiter="\t")
            for line_number, row in enumerate(raw_reader, start=1):
                if not row:
                    continue
                if row[0].startswith("#"):
                    continue
                if len(row) < 4:
                    raise ValueError(
                        f"Positions file '{path}' has fewer than 4 columns on line {line_number}."
                    )

                key = row[0].strip()
                chrom = row[1].strip()
                start = row[2].strip()
                end = row[3].strip()
                species = row[4].strip() if len(row) >= 5 else ""

                if not key or not chrom or not start or not end:
                    raise ValueError(
                        f"Positions file '{path}' has an empty required field on line {line_number}."
                    )

                new_value = {
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "species": species,
                }

                if key in pos_dict:
                    if pos_dict[key] != new_value:
                        raise ValueError(
                            f"Conflicting position information for key '{key}' in positions "
                            f"file '{path}' on line {line_number}. First value: {pos_dict[key]}; "
                            f"new value: {new_value}"
                        )
                else:
                    pos_dict[key] = new_value

    if not pos_dict:
        raise ValueError(f"No positions were loaded from: {path}")

    return pos_dict


# ----------------------------------------------------------------------
# Classification logic
# ----------------------------------------------------------------------

def classify_species_set(
    species_set: set[str],
    branch_to_species: Dict[str, set[str]]
) -> str:
    """
    Return the smallest branch_id whose species set contains all species in
    species_set.
    """
    for branch_id in sorted(branch_to_species.keys(), key=branch_sort_key):
        if species_set.issubset(branch_to_species[branch_id]):
            return branch_id

    raise ValueError(
        "Species set did not match any branch definition: "
        f"{','.join(sorted(species_set))}"
    )


def build_output_row(
    row_dict: Dict[str, str],
    pos_dict: Dict[str, Dict[str, str]],
    branch_to_species: Dict[str, set[str]],
    color_dict: Dict[str, str],
    tip_separator: str,
    label_format: str,
    position_key_type: str,
    include_metadata: bool
) -> List[str]:
    target_tips_field = row_dict["target_tips"].strip()
    sharing_species_field = row_dict["sharing_species"].strip()

    target_tips = [x.strip() for x in target_tips_field.split(",") if x.strip()]
    if len(target_tips) != 2:
        raise ValueError(
            f"Expected exactly 2 target tips, found {len(target_tips)}: '{target_tips_field}'"
        )

    parsed_tip1 = parse_tip_label(target_tips[0], tip_separator, label_format)
    parsed_tip2 = parse_tip_label(target_tips[1], tip_separator, label_format)

    tip_species1 = parsed_tip1["species"]
    tip_species2 = parsed_tip2["species"]

    sharing_species = {x.strip() for x in sharing_species_field.split(",") if x.strip()}
    if not sharing_species:
        raise ValueError("sharing_species field is empty.")

    redip_branch = classify_species_set(sharing_species, branch_to_species)
    color = color_dict[redip_branch]

    key1 = extract_position_key(
        tip_label=target_tips[0],
        tip_separator=tip_separator,
        label_format=label_format,
        position_key_type=position_key_type
    )
    key2 = extract_position_key(
        tip_label=target_tips[1],
        tip_separator=tip_separator,
        label_format=label_format,
        position_key_type=position_key_type
    )

    if key1 not in pos_dict:
        raise ValueError(f"Position key '{key1}' was not found in the positions table.")
    if key2 not in pos_dict:
        raise ValueError(f"Position key '{key2}' was not found in the positions table.")

    pos_species1 = pos_dict[key1].get("species", "")
    pos_species2 = pos_dict[key2].get("species", "")

    if pos_species1 and pos_species1 != tip_species1:
        raise ValueError(
            f"Species mismatch for key '{key1}': positions file has "
            f"'{pos_species1}', tip label has '{tip_species1}'."
        )

    if pos_species2 and pos_species2 != tip_species2:
        raise ValueError(
            f"Species mismatch for key '{key2}': positions file has "
            f"'{pos_species2}', tip label has '{tip_species2}'."
        )

    chr1 = pos_dict[key1]["chr"]
    start1 = pos_dict[key1]["start"]
    end1 = pos_dict[key1]["end"]

    chr2 = pos_dict[key2]["chr"]
    start2 = pos_dict[key2]["start"]
    end2 = pos_dict[key2]["end"]

    output_row = [chr1, start1, end1, chr2, start2, end2, color]

    if include_metadata:
        output_row.extend([
            redip_branch,
            row_dict.get("tree_file", ""),
            row_dict.get("tree_index", ""),
            row_dict.get("allowed_clade_index", ""),
            row_dict.get("num_tips_in_allowed_clade", ""),
            target_tips_field,
            sharing_species_field,
        ])

    return output_row


# ----------------------------------------------------------------------
# File I/O
# ----------------------------------------------------------------------

def prepare_output_file(path: Path, on_exists: str) -> str:
    if path.exists():
        if on_exists == "error":
            raise FileExistsError(f"Output file already exists: {path}")
        if on_exists == "overwrite":
            return "w"
        if on_exists == "append":
            return "a"

    return "w"


def get_input_rows(
    classification_path: str,
    no_header: bool
) -> Iterable[Dict[str, str]]:
    """
    Yield classification rows as dictionaries.

    Supports:
    - headered input from the rewritten script 1
    - no-header input in the same column order
    """
    expected_columns = [
        "tree_file",
        "tree_index",
        "allowed_clade_index",
        "num_tips_in_allowed_clade",
        "target_tips",
        "sharing_species",
    ]

    with open(classification_path, "r", newline="") as handle:
        if no_header:
            reader = csv.reader(handle, delimiter="\t")
            for line_number, row in enumerate(reader, start=1):
                if not row:
                    continue
                if len(row) < 6:
                    raise ValueError(
                        f"Classification file '{classification_path}' has fewer than 6 "
                        f"columns on line {line_number}."
                    )
                yield dict(zip(expected_columns, row[:6]))
        else:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(
                    f"Classification file '{classification_path}' appears to have no header."
                )

            missing = [col for col in expected_columns if col not in reader.fieldnames]
            if missing:
                raise ValueError(
                    f"Classification file '{classification_path}' is missing required "
                    f"columns: {', '.join(missing)}"
                )

            for row in reader:
                yield {col: (row[col] or "") for col in expected_columns}


def write_rows(
    output_path: str,
    rows: Iterable[List[str]],
    write_header: bool,
    include_metadata: bool,
    mode: str
) -> int:
    rows = list(rows)

    header = ["chr1", "start1", "end1", "chr2", "start2", "end2", "color"]
    if include_metadata:
        header.extend([
            "redip_branch",
            "tree_file",
            "tree_index",
            "allowed_clade_index",
            "num_tips_in_allowed_clade",
            "target_tips",
            "sharing_species",
        ])

    should_write_header = False
    if write_header:
        if mode == "w":
            should_write_header = True
        elif mode == "a":
            output_file = Path(output_path)
            if (not output_file.exists()) or output_file.stat().st_size == 0:
                should_write_header = True

    with open(output_path, mode, newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        if should_write_header:
            writer.writerow(header)

        for row in rows:
            writer.writerow(row)

    return len(rows)


def process_classification_file(args: argparse.Namespace) -> List[List[str]]:
    branch_to_species = load_branch_definitions(args.branch_definitions)
    color_dict = load_colors(args.colors, list(branch_to_species.keys()))
    pos_dict = load_positions(
        path=args.positions,
        has_header=args.position_has_header,
        key_column=args.position_key_column,
        chr_column=args.position_chr_column,
        start_column=args.position_start_column,
        end_column=args.position_end_column,
        species_column=args.position_species_column
    )

    output_rows: List[List[str]] = []

    for record_index, row_dict in enumerate(
        get_input_rows(args.classification_tsv, args.classification_no_header),
        start=1
    ):
        try:
            output_row = build_output_row(
                row_dict=row_dict,
                pos_dict=pos_dict,
                branch_to_species=branch_to_species,
                color_dict=color_dict,
                tip_separator=args.tip_separator,
                label_format=args.label_format,
                position_key_type=args.position_key_type,
                include_metadata=args.include_metadata
            )
            output_rows.append(output_row)

        except Exception as exc:
            message = (
                f"Failed while processing classification record {record_index} "
                f"from {args.classification_tsv}: {exc}"
            )
            handle_error(message, args.strict)
            continue

    return output_rows


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    output_path = Path(args.output)
    mode = prepare_output_file(output_path, args.on_exists)

    rows = process_classification_file(args)

    write_rows(
        output_path=args.output,
        rows=rows,
        write_header=args.write_header,
        include_metadata=args.include_metadata,
        mode=mode
    )

    logging.info("Finished.")
    logging.info("Output file: %s", args.output)
    logging.info("Rows written: %d", len(rows))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.error(str(exc))
        sys.exit(1)
