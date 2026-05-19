#!/usr/bin/env python3

"""
Shared helpers for the rediploidisation workflow.

This module centralises parsing of genomes.tsv so all redip scripts interpret
species IDs, focal-species inclusion, per-species classification modes, and
outgroup flags consistently.

Backward compatibility notes
----------------------------
Older genomes.tsv files used a rediploidisation column with 1/0 values.
That remains supported:
    1 -> run redip using the pipeline default mode
    0 -> do not run redip for that species

Newer genomes.tsv files can instead use text values in the same column:
    standard
    recurrent
    recurrent_lossy
    none
"""

from __future__ import annotations
import csv
from pathlib import Path
from typing import Iterable, List, Tuple


SPECIES_HEADER_NAMES = {"species", "species_id", "taxon", "genome_id", "genome"}
REDIP_HEADER_NAMES = {
    "rediploidisation",
    "rediploidization",
    "redip",
    "wgd",
    "include_redip",
    "redip_ingroup",
    "redip_mode",
}
OUTGROUP_HEADER_NAMES = {"outgroup", "is_outgroup", "root_outgroup"}

# Kept for backward compatibility with older helper functions that only needed
# a boolean redip inclusion decision.
TRUE_VALUES = {
    "true",
    "t",
    "yes",
    "y",
    "1",
    "wgd",
    "redip",
    "rediploidisation",
    "rediploidization",
}

FALSE_VALUES = {"", "false", "f", "no", "n", "0", "none", "skip", "na", "nan"}

# Accepted per-species redip mode values in genomes.tsv.
REDIP_SKIP_VALUES = FALSE_VALUES
REDIP_DEFAULT_VALUES = TRUE_VALUES | {"default", "run"}
REDIP_STANDARD_VALUES = {"standard", "exact", "exact_n", "standard_exact_n"}
REDIP_RECURRENT_VALUES = {"recurrent", "balanced"}
REDIP_LOSSY_VALUES = {"recurrent_lossy", "recurrent-lossy", "lossy"}
VALID_GLOBAL_REDIP_MODES = {"standard", "recurrent", "recurrent_lossy"}
VALID_NORMALISED_REDIP_MODES = {"none", "default", *VALID_GLOBAL_REDIP_MODES}


def dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _find_column(fieldnames: list[str], accepted: set[str]) -> str | None:
    """Return the actual column name matching any accepted lowercase alias."""
    normalized = {name.strip().lower(): name for name in fieldnames if name is not None}
    for name in accepted:
        if name in normalized:
            return normalized[name]
    return None


def normalise_redip_mode(value: object, default_mode: str = "recurrent") -> str:
    """
    Convert a genomes.tsv redip value to one normalised analysis mode.

    The 1/0 behaviour is intentionally preserved for backward compatibility:
      - 1, true, yes, redip, wgd -> default_mode
      - 0, false, no, none       -> none

    Text values allow per-species overrides:
      - standard
      - recurrent
      - recurrent_lossy
    """
    default_mode = str(default_mode).strip().lower()
    if default_mode not in VALID_GLOBAL_REDIP_MODES:
        raise ValueError(
            "Invalid default redip classification mode: "
            f"{default_mode}. Expected one of: "
            + ", ".join(sorted(VALID_GLOBAL_REDIP_MODES))
        )

    mode = "" if value is None else str(value).strip().lower()

    if mode in REDIP_SKIP_VALUES:
        return "none"
    if mode in REDIP_DEFAULT_VALUES:
        return default_mode
    if mode in REDIP_STANDARD_VALUES:
        return "standard"
    if mode in REDIP_RECURRENT_VALUES:
        return "recurrent"
    if mode in REDIP_LOSSY_VALUES:
        return "recurrent_lossy"

    raise ValueError(
        f"Unknown rediploidisation value: '{value}'. "
        "Use 0/1 for backward-compatible behaviour or one of: "
        "none, default, standard, recurrent, recurrent_lossy."
    )


def read_species_list(path: str) -> List[str]:
    """
    Read a standalone allowed-species TSV.

    Recognises headers:
        species, species_id, taxon, genome_id, genome

    If no recognised header exists, uses the first column. This is kept for
    backward compatibility with older --allowed-species files.
    """
    rows = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row and not row[0].startswith("#")]

    if not rows:
        raise ValueError(f"Species file is empty: {path}")

    first = [x.strip().lower() for x in rows[0]]
    species_col_idx = None

    for i, value in enumerate(first):
        if value in SPECIES_HEADER_NAMES:
            species_col_idx = i
            break

    data_rows = rows[1:] if species_col_idx is not None else rows
    if species_col_idx is None:
        species_col_idx = 0

    species = []
    for row in data_rows:
        if species_col_idx >= len(row):
            continue
        value = row[species_col_idx].strip()
        if value:
            species.append(value)

    species = dedupe_preserve_order(species)

    if not species:
        raise ValueError(f"No species were loaded from: {path}")

    return species


def _detect_genomes_tsv_columns(path: str) -> tuple[list[str], str, str]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"genomes.tsv appears to have no header: {path}")

        species_col = _find_column(reader.fieldnames, SPECIES_HEADER_NAMES)
        redip_col = _find_column(reader.fieldnames, REDIP_HEADER_NAMES)

    if species_col is None:
        raise ValueError(
            "Could not find a species/genome column in genomes.tsv. "
            f"Accepted names: {', '.join(sorted(SPECIES_HEADER_NAMES))}"
        )

    if redip_col is None:
        raise ValueError(
            "Could not find a rediploidisation/WGD inclusion column in genomes.tsv. "
            f"Accepted names: {', '.join(sorted(REDIP_HEADER_NAMES))}"
        )

    return reader.fieldnames or [], species_col, redip_col


def read_redip_species_modes_from_genomes_tsv(
    path: str,
    default_mode: str = "recurrent",
) -> list[tuple[str, str]]:
    """
    Return focal redip species and their resolved classification mode.

    Output rows are (species_id, mode), where mode is one of:
        standard, recurrent, recurrent_lossy

    Rows resolving to mode 'none' are skipped.
    """
    _, species_col, redip_col = _detect_genomes_tsv_columns(path)

    rows: list[tuple[str, str]] = []
    seen = set()

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for line_number, row in enumerate(reader, start=2):
            species_id = (row.get(species_col) or "").strip()
            raw_value = row.get(redip_col)

            if not species_id:
                raise ValueError(
                    f"Missing species/genome value on line {line_number} of {path}"
                )

            mode = normalise_redip_mode(raw_value, default_mode=default_mode)
            if mode == "none":
                continue

            if species_id not in seen:
                seen.add(species_id)
                rows.append((species_id, mode))

    if not rows:
        raise ValueError(
            f"No rediploidisation species were selected from genomes.tsv: {path}"
        )

    return rows


def read_redip_species_from_genomes_tsv(
    path: str,
    default_mode: str = "recurrent",
) -> List[str]:
    """
    Read rediploidisation/WGD ingroup species from genomes.tsv.

    Kept for backward compatibility with existing scripts. This now treats any
    non-'none' resolved redip mode as included and returns species IDs only.
    """
    return [
        species
        for species, _mode in read_redip_species_modes_from_genomes_tsv(
            path,
            default_mode=default_mode,
        )
    ]


def read_redip_species(
    genomes_tsv: str | None = None,
    allowed_species: str | None = None,
    default_mode: str = "recurrent",
) -> List[str]:
    """
    Prefer genomes.tsv if supplied. Fall back to standalone allowed-species file.

    The default_mode argument is ignored for standalone species lists because
    those files contain only allowed species, not per-species mode values.
    """
    if genomes_tsv:
        return read_redip_species_from_genomes_tsv(
            genomes_tsv,
            default_mode=default_mode,
        )

    if allowed_species:
        return read_species_list(allowed_species)

    raise ValueError("Provide either --genomes-tsv or --allowed-species.")


def read_outgroup_species_from_genomes_tsv(path: str) -> list[str]:
    """
    Read outgroup species from genomes.tsv.

    This preserves the earlier boolean-style behaviour. Tiered outgroup rooting
    is handled by root_gene_tree.py, which reads the outgroup column directly.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"genomes.tsv appears to have no header: {path}")

        species_col = _find_column(reader.fieldnames, SPECIES_HEADER_NAMES)
        outgroup_col = _find_column(reader.fieldnames, OUTGROUP_HEADER_NAMES)

        if species_col is None:
            raise ValueError(
                "Could not find a species/genome column in genomes.tsv. "
                f"Accepted names: {', '.join(sorted(SPECIES_HEADER_NAMES))}"
            )

        if outgroup_col is None:
            raise ValueError(
                "Could not find an outgroup column in genomes.tsv. "
                f"Accepted names: {', '.join(sorted(OUTGROUP_HEADER_NAMES))}"
            )

        species = []
        for line_number, row in enumerate(reader, start=2):
            species_id = (row.get(species_col) or "").strip()
            include_value = (row.get(outgroup_col) or "").strip().lower()

            if not species_id:
                raise ValueError(
                    f"Missing species/genome value on line {line_number} of {path}"
                )

            if include_value in TRUE_VALUES or include_value == "outgroup":
                species.append(species_id)

    species = dedupe_preserve_order(species)

    if not species:
        raise ValueError(f"No outgroup species were selected from genomes.tsv: {path}")

    return species
