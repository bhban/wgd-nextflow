#!/usr/bin/env python3

from __future__ import annotations

"""
Prepare Circos input directories for one species and one plot layer.

Input is the metadata-rich links TSV from make_links.py. This script writes:

    complete/
        karyotype.txt
        <plot>.links.tsv
        circos.conf

    branch_specific/branch_<redip_branch>/
        karyotype.txt
        <plot>.branch_<redip_branch>.links.tsv
        circos.conf

The complete plot shows all links for the species/layer. Branch-specific plots
show only links assigned to one redip branch/colour, which helps with visual
clutter. This replaces the older behaviour where per-branch link files were
written but only one top-level circos.conf was generated.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare complete and branch-specific Circos inputs."
    )
    parser.add_argument("--species", required=True, help="Plot label, usually species.plot_level.")
    parser.add_argument("--circos-links", required=True)
    parser.add_argument("--chr-bed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--allow-empty",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "If no links remain after filtering, write a README and branch_summary.tsv "
            "but do not create circos.conf files. Enabled by default so empty plot "
            "levels do not fail the whole workflow."
        ),
    )
    return parser.parse_args()


def safe_label(value: str) -> str:
    """Make a value safe for filenames while preserving readable branch IDs."""
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return value or "unknown"


def branch_sort_key(value: str) -> tuple[int, str]:
    try:
        return (0, f"{int(value):010d}")
    except ValueError:
        return (1, value)



def chromosome_sort_key(chrom: str):
    """
    Sort chromosomes naturally, grouping subgenomes.

    Examples:
      chr1, chr2, chr10
      chrA1, chrA2, chrA10, chrB1, chrB2, chrB10
    """
    match = re.fullmatch(r"chr([A-Za-z]*)(\d+)", chrom, re.IGNORECASE)

    if match:
        subgenome = match.group(1).upper()
        number = int(match.group(2))
        return (0, subgenome, number, chrom)

    parts = re.split(r"(\d+)", chrom)
    natural = tuple(
        (0, int(part)) if part.isdigit()
        else (1, part.lower())
        for part in parts
    )
    return (1, natural)


def make_karyotype(chr_bed: Path, out_karyotype: Path) -> list[str]:
    """
    Create a Circos karyotype from the chromosome-length BED.
    Preserve original chromosome IDs and sort naturally.
    """
    chromosomes = {}
    out_karyotype.parent.mkdir(parents=True, exist_ok=True)

    with open(chr_bed, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 3:
                continue

            chrom = row[0].strip()

            try:
                start = int(row[1])
                end = int(row[2])
            except ValueError:
                raise ValueError(f"Invalid chromosome coordinates: {row}")

            if start < 0 or end <= start:
                raise ValueError(f"Invalid chromosome interval: {row}")

            if chrom in chromosomes:
                raise ValueError(f"Duplicate chromosome in BED: {chrom}")

            chromosomes[chrom] = (start, end)

    if not chromosomes:
        raise ValueError(f"No valid chromosomes found in {chr_bed}")

    chr_names = sorted(chromosomes, key=chromosome_sort_key)

    with open(out_karyotype, "w") as out:
        for chrom in chr_names:
            start, end = chromosomes[chrom]
            out.write(
                f"chr\t-\t{chrom}\t{chrom}\t{start}\t{end}\tblack\n"
            )

    return chr_names



def load_links(circos_tsv: Path) -> list[dict[str, str]]:
    with open(circos_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Circos links file appears to have no header: {circos_tsv}")

        required = {
            "chr1", "start1", "end1",
            "chr2", "start2", "end2",
            "color", "redip_branch",
        }

        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(
                f"Circos links file {circos_tsv} is missing required columns: "
                + ", ".join(missing)
            )

        return list(reader)


def filter_main_chr_rows(rows: list[dict[str, str]], valid_chr: set[str]) -> list[dict[str, str]]:
    valid_chr_lookup = {chrom.lower(): chrom for chrom in valid_chr}
    kept: list[dict[str, str]] = []
    for row in rows:
        chr1 = (row.get("chr1") or "").strip()
        chr2 = (row.get("chr2") or "").strip()
        chr1_key = chr1.lower()
        chr2_key = chr2.lower()
        if chr1_key in valid_chr_lookup and chr2_key in valid_chr_lookup:
            updated_row = row.copy()
            updated_row["chr1"] = valid_chr_lookup[chr1_key]
            updated_row["chr2"] = valid_chr_lookup[chr2_key]
            kept.append(updated_row)
    return kept


def write_links_file(rows: list[dict[str, str]], out_path: Path) -> None:
    """Write the 7-column Circos link file consumed by circos.conf."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        for row in rows:
            writer.writerow([
                row["chr1"].strip(),
                row["start1"].strip(),
                row["end1"].strip(),
                row["chr2"].strip(),
                row["start2"].strip(),
                row["end2"].strip(),
                row["color"].strip(),
            ])


def write_circos_conf(
    plot_name: str,
    out_conf: Path,
    karyotype_path: Path,
    links_path: Path,
) -> None:
    """Write a self-contained Circos configuration for one links file."""
    conf = f"""karyotype = {karyotype_path.name}

chromosomes_units = 1000000

<ideogram>
<spacing>
default = 0.01r
</spacing>

radius           = 0.85r
thickness        = 30p
fill             = yes
stroke_color     = black
stroke_thickness = 2p

show_label       = yes
label_font       = default
label_radius     = dims(ideogram,radius_outer) + 0.08r
label_size       = 64p
label_parallel   = yes
</ideogram>

<ticks>
show_ticks       = yes
show_tick_labels = yes

<tick>
spacing        = 10u
size           = 8p
thickness      = 2p
color          = black
show_label     = no
</tick>

<tick>
spacing        = 50u
size           = 12p
thickness      = 2p
color          = black
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>
</ticks>

<links>
<link>
file          = {links_path.name}
radius        = 0.78r
bezier_radius = 0.1r
thickness     = 2p
color         = eval(var(color))
z             = 10
record_limit  = 100000
</link>
</links>

<image>
dir               = .
file              = circos_{plot_name}
png               = yes
svg               = yes
radius            = 1800p
background        = white
angle_offset      = -90
auto_alpha_colors = yes
auto_alpha_steps  = 5
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
"""
    out_conf.write_text(conf)


def prepare_one_circos_dir(
    plot_name: str,
    out_dir: Path,
    chr_bed: Path,
    rows: list[dict[str, str]],
) -> None:
    """Write one runnable Circos directory."""
    out_dir.mkdir(parents=True, exist_ok=True)
    karyotype_path = out_dir / "karyotype.txt"
    links_path = out_dir / f"{plot_name}.links.tsv"
    conf_path = out_dir / "circos.conf"

    make_karyotype(chr_bed, karyotype_path)
    write_links_file(rows, links_path)
    write_circos_conf(
        plot_name=plot_name,
        out_conf=conf_path,
        karyotype_path=karyotype_path,
        links_path=links_path,
    )

    with open(out_dir / "README.txt", "w") as out:
        out.write(f"Plot name: {plot_name}\n")
        out.write(f"Rows: {len(rows)}\n")
        out.write(f"Links: {links_path.name}\n")
        out.write(f"Circos config: {conf_path.name}\n")


def write_empty_outputs(output_dir: Path, species: str, reason: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "README.txt", "w") as out:
        out.write(f"Species/plot: {species}\n")
        out.write("No Circos config was written.\n")
        out.write(f"Reason: {reason}\n")
    with open(output_dir / "branch_summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["branch_id", "link_rows"])
        writer.writerow(["complete", 0])
    (output_dir / "NO_LINKS.txt").write_text(reason + "\n")


def main() -> None:
    args = parse_args()

    species = args.species
    circos_tsv = Path(args.circos_links)
    chr_bed = Path(args.chr_bed)
    output_dir = Path(args.output_dir)

    if not chr_bed.exists():
        raise FileNotFoundError(f"Missing chromosome BED for {species}: {chr_bed}")
    if not circos_tsv.exists():
        raise FileNotFoundError(f"Missing Circos links TSV for {species}: {circos_tsv}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build a temporary karyotype once to identify valid main chromosomes.
    tmp_karyotype = output_dir / "_tmp_karyotype_check.txt"
    chr_names = make_karyotype(chr_bed, tmp_karyotype)
    tmp_karyotype.unlink(missing_ok=True)
    valid_chr = set(chr_names)

    all_rows = load_links(circos_tsv)
    kept_rows = filter_main_chr_rows(all_rows, valid_chr)

    if not kept_rows:
        reason = f"No links remained after filtering to main chromosomes for {species}."
        if args.allow_empty:
            write_empty_outputs(output_dir, species, reason)
            print(reason)
            return
        raise ValueError(reason)

    # Complete plot: all links for this species/layer.
    complete_dir = output_dir / "complete"
    prepare_one_circos_dir(
        plot_name=safe_label(species),
        out_dir=complete_dir,
        chr_bed=chr_bed,
        rows=kept_rows,
    )

    # Branch-specific plots: one plot per redip_branch/colour.
    # Empty branch IDs are treated as an error because those links would appear
    # in the complete plot but silently disappear from branch-specific plots.
    branch_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row_number, row in enumerate(kept_rows, start=1):
        branch = (row.get("redip_branch") or "").strip()
        if not branch:
            raise ValueError(
                f"Row {row_number} in {circos_tsv} has an empty redip_branch. "
                "Cannot create branch-specific Circos plots."
            )
        branch_rows[branch].append(row)

    branch_specific_root = output_dir / "branch_specific"
    for branch in sorted(branch_rows, key=branch_sort_key):
        branch_label = f"branch_{safe_label(branch)}"
        plot_name = f"{safe_label(species)}.{branch_label}"
        prepare_one_circos_dir(
            plot_name=plot_name,
            out_dir=branch_specific_root / branch_label,
            chr_bed=chr_bed,
            rows=branch_rows[branch],
        )

    with open(output_dir / "branch_summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["branch_id", "link_rows"])
        writer.writerow(["complete", len(kept_rows)])
        for branch in sorted(branch_rows, key=branch_sort_key):
            writer.writerow([branch, len(branch_rows[branch])])

    with open(output_dir / "README.txt", "w") as out:
        out.write(f"Species/plot: {species}\n")
        out.write(f"Rows retained on main chromosomes: {len(kept_rows)}\n")
        out.write("Complete plot directory: complete\n")
        out.write("Branch-specific plot directories:\n")
        for branch in sorted(branch_rows, key=branch_sort_key):
            branch_label = f"branch_{safe_label(branch)}"
            out.write(
                f"  branch_specific/{branch_label}\t"
                f"{len(branch_rows[branch])} rows\n"
            )

    print(f"Prepared Circos inputs for {species}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
