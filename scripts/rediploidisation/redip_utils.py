from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List


SPECIES_HEADER_NAMES = {"species", "species_id", "taxon", "genome_id", "genome"}
REDIP_HEADER_NAMES = {
    "rediploidisation",
    "rediploidization",
    "redip",
    "wgd",
    "include_redip",
    "redip_ingroup",
}

OUTGROUP_HEADER_NAMES = {"outgroup", "is_outgroup", "root_outgroup"}

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


def dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _find_column(fieldnames: list[str], accepted: set[str]) -> str | None:
    normalized = {name.strip().lower(): name for name in fieldnames if name is not None}
    for name in accepted:
        if name in normalized:
            return normalized[name]
    return None


def read_species_list(path: str) -> List[str]:
    """
    Read a standalone allowed-species TSV.

    Recognises headers:
        species, species_id, taxon, genome_id, genome

    If no recognised header exists, uses the first column.
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


def read_redip_species_from_genomes_tsv(path: str) -> List[str]:
    """
    Read rediploidisation/WGD ingroup species from genomes.tsv.

    Requires:
      - a genome/species column recognised as one of:
        genome, genome_id, species, species_id, taxon
      - an optional rediploidisation column recognised as one of:
        rediploidisation, rediploidization, redip, wgd,
        include_redip, redip_ingroup

    Rows where the redip column is one of TRUE_VALUES are included.
    """
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

        species = []
        for line_number, row in enumerate(reader, start=2):
            species_id = (row.get(species_col) or "").strip()
            include_value = (row.get(redip_col) or "").strip().lower()

            if not species_id:
                raise ValueError(
                    f"Missing species/genome value on line {line_number} of {path}"
                )

            if include_value in TRUE_VALUES:
                species.append(species_id)

    species = dedupe_preserve_order(species)

    if not species:
        raise ValueError(
            f"No rediploidisation species were selected from genomes.tsv: {path}"
        )

    return species


def read_redip_species(
    genomes_tsv: str | None = None,
    allowed_species: str | None = None,
) -> List[str]:
    """
    Prefer genomes.tsv if supplied. Fall back to standalone allowed-species file.
    """
    if genomes_tsv:
        return read_redip_species_from_genomes_tsv(genomes_tsv)

    if allowed_species:
        return read_species_list(allowed_species)

    raise ValueError(
        "Provide either --genomes-tsv or --allowed-species."
    )


def read_outgroup_species_from_genomes_tsv(path: str) -> list[str]:
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
