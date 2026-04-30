#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from ete3 import Tree


FALSE_VALUES = {"", "false", "f", "no", "n", "0", "none", "na", "nan"}
TRUE_VALUES = {"true", "t", "yes", "y"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Root a gene tree on outgroup species specified in genomes.tsv. "
            "The outgroup column may be boolean or integer-tiered. "
            "For integer tiers, lower values are treated as more basal and tried first."
        )
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree.")
    parser.add_argument("--genomes-tsv", required=True, help="Pipeline genomes.tsv.")
    parser.add_argument("--output-tree", required=True, help="Output rooted Newick tree.")
    parser.add_argument("--summary-tsv", required=True, help="Rooting summary TSV.")
    parser.add_argument("--tip-separator", default="|", help="Tip-label separator.")
    parser.add_argument("--tree-format", type=int, default=1, help="ETE tree format.")
    parser.add_argument(
        "--species-column",
        default=None,
        help=(
            "Column in genomes.tsv containing species IDs. "
            "If omitted, tries species, species_id, taxon, then genome_id."
        ),
    )
    parser.add_argument(
        "--outgroup-column",
        default="outgroup",
        help=(
            "Column in genomes.tsv containing outgroup status or tier. "
            "Values may be TRUE/FALSE or integer tiers."
        ),
    )
    return parser.parse_args()


def get_species(tip_name: str, separator: str) -> str:
    return tip_name.split(separator, 1)[0]


def parse_outgroup_value(raw_value: object) -> int | None:
    """
    Convert an outgroup value from genomes.tsv into an outgroup tier.

    Returns
    -------
    int | None
        Integer tier if the row is an outgroup.
        None if the row is not an outgroup.

    Accepted examples
    -----------------
    TRUE, true, yes, y -> tier 1
    FALSE, false, no, n, 0, blank -> not outgroup
    1, 2, 3 -> corresponding tier
    """

    value = "" if raw_value is None else str(raw_value).strip()
    value_lower = value.lower()

    if value_lower in FALSE_VALUES:
        return None

    if value_lower in TRUE_VALUES:
        return 1

    try:
        tier = int(value)
    except ValueError as exc:
        raise ValueError(
            f"Invalid outgroup value '{value}'. "
            "Use TRUE/FALSE or positive integer tiers."
        ) from exc

    if tier <= 0:
        return None

    return tier


def detect_species_column(
    fieldnames: list[str],
    requested_species_column: str | None,
) -> str:
    if requested_species_column:
        if requested_species_column not in fieldnames:
            raise ValueError(
                f"Requested species column '{requested_species_column}' was not found. "
                f"Available columns: {', '.join(fieldnames)}"
            )
        return requested_species_column

    candidates = ["species", "species_id", "taxon", "genome_id"]

    for candidate in candidates:
        if candidate in fieldnames:
            return candidate

    raise ValueError(
        "Could not detect species column in genomes.tsv. "
        "Expected one of: species, species_id, taxon, genome_id. "
        "Alternatively, provide --species-column."
    )


def read_outgroup_tiers_from_genomes_tsv(
    genomes_tsv: str,
    outgroup_column: str,
    species_column: str | None = None,
) -> dict[int, set[str]]:
    """
    Read outgroup species from genomes.tsv as tiered groups.

    Boolean TRUE values are treated as tier 1, preserving the previous behaviour
    when there is only one outgroup clade.
    """

    tiers: dict[int, set[str]] = defaultdict(set)

    with open(genomes_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Could not read header from {genomes_tsv}")

        fieldnames = reader.fieldnames

        if outgroup_column not in fieldnames:
            raise ValueError(
                f"Outgroup column '{outgroup_column}' was not found in {genomes_tsv}. "
                f"Available columns: {', '.join(fieldnames)}"
            )

        detected_species_column = detect_species_column(
            fieldnames=fieldnames,
            requested_species_column=species_column,
        )

        for row_number, row in enumerate(reader, start=2):
            species = str(row.get(detected_species_column, "")).strip()

            if not species:
                raise ValueError(
                    f"Missing species value in column '{detected_species_column}' "
                    f"on line {row_number}"
                )

            tier = parse_outgroup_value(row.get(outgroup_column))

            if tier is not None:
                tiers[tier].add(species)

    return dict(tiers)


def root_on_candidate_species(
    tree: Tree,
    candidate_species: set[str],
    tip_separator: str,
) -> tuple[bool, str]:
    """
    Try rooting on one outgroup tier.

    If one matching tip is present, root on that tip.
    If multiple matching tips are present, root on their MRCA.
    If their MRCA is already the whole tree, leave the tree unchanged.
    """

    candidate_leaves = [
        leaf for leaf in tree.iter_leaves()
        if get_species(leaf.name, tip_separator) in candidate_species
    ]

    n_candidate_tips = len(candidate_leaves)

    if n_candidate_tips == 0:
        return False, "no_tips_present"

    found_species = sorted({
        get_species(leaf.name, tip_separator)
        for leaf in candidate_leaves
    })

    if n_candidate_tips == 1:
        tree.set_outgroup(candidate_leaves[0])
        return True, (
            "rooted_on_single_outgroup_tip"
            + "__species="
            + found_species[0]
        )

    mrca = tree.get_common_ancestor(candidate_leaves)

    if mrca is tree:
        return False, (
            "outgroup_mrca_is_root__tree_left_unchanged"
            + "__n_outgroup_tips="
            + str(n_candidate_tips)
            + "__species="
            + ",".join(found_species)
        )

    tree.set_outgroup(mrca)

    return True, (
        "rooted_on_outgroup_mrca"
        + "__n_outgroup_tips="
        + str(n_candidate_tips)
        + "__species="
        + ",".join(found_species)
    )


def root_on_tiered_outgroups(
    tree: Tree,
    outgroup_tiers: dict[int, set[str]],
    tip_separator: str,
) -> tuple[bool, str]:
    """
    Try outgroup tiers from most basal to least basal.

    Lower tier numbers are considered more basal.

    For each tier:
      - skip if no species from that tier are present in the gene tree
      - root on a single present tip if only one matching tip exists
      - otherwise root on the MRCA of all matching tips
      - if the MRCA is the whole tree, try the next tier inward
    """

    if not outgroup_tiers:
        return False, "no_outgroup_species_defined"

    attempt_messages: list[str] = []

    for tier in sorted(outgroup_tiers):
        candidate_species = outgroup_tiers[tier]

        rooted, message = root_on_candidate_species(
            tree=tree,
            candidate_species=candidate_species,
            tip_separator=tip_separator,
        )

        tier_message = (
            "tier="
            + str(tier)
            + "__defined_species="
            + ",".join(sorted(candidate_species))
            + "__"
            + message
        )

        attempt_messages.append(tier_message)

        if rooted:
            return True, tier_message

    return False, "no_rootable_outgroup_tier__" + " || ".join(attempt_messages)


def write_summary(path: str, tree_id: str, status: str, message: str) -> None:
    clean_message = message.replace("\t", " ").replace("\n", " ")

    with open(path, "w") as out:
        out.write("tree_file\tstatus\tmessage\n")
        out.write(f"{tree_id}\t{status}\t{clean_message}\n")


def main() -> None:
    args = parse_args()

    tree_id = Path(args.tree).name

    try:
        outgroup_tiers = read_outgroup_tiers_from_genomes_tsv(
            genomes_tsv=args.genomes_tsv,
            outgroup_column=args.outgroup_column,
            species_column=args.species_column,
        )
    except Exception as exc:
        write_summary(
            args.summary_tsv,
            tree_id,
            "ERROR",
            f"failed_to_read_outgroups__{exc}",
        )
        sys.exit(1)

    try:
        tree = Tree(args.tree, format=args.tree_format)
    except Exception as exc:
        write_summary(
            args.summary_tsv,
            tree_id,
            "ERROR",
            f"failed_to_read_tree__{exc}",
        )
        sys.exit(1)

    try:
        rooted, message = root_on_tiered_outgroups(
            tree=tree,
            outgroup_tiers=outgroup_tiers,
            tip_separator=args.tip_separator,
        )

        tree.write(outfile=args.output_tree, format=args.tree_format)

        status = "ROOTED" if rooted else "SKIPPED"
        write_summary(args.summary_tsv, tree_id, status, message)

    except Exception as exc:
        write_summary(
            args.summary_tsv,
            tree_id,
            "ERROR",
            f"rooting_failed__{exc}",
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
