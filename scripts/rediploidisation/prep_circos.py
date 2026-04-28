#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Circos input files for one rediploidisation species."
    )
    parser.add_argument("--species", required=True)
    parser.add_argument("--circos-links", required=True)
    parser.add_argument("--chr-bed", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def branch_sort_key(value: str) -> tuple[int, str]:
    try:
        return (0, f"{int(value):010d}")
    except ValueError:
        return (1, value)


def make_karyotype(chr_bed: Path, out_karyotype: Path) -> list[str]:
    chr_names = []

    with open(chr_bed, newline="") as handle, open(out_karyotype, "w", newline="") as out:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row:
                continue

            if row[0].startswith("#"):
                continue

            if len(row) < 3:
                continue

            chrom = row[0].strip()
            start = row[1].strip()
            end = row[2].strip()

            if not re.fullmatch(r"chr[0-9]+", chrom):
                continue

            out.write(f"chr\t-\t{chrom}\t{chrom}\t{start}\t{end}\tblack\n")
            chr_names.append(chrom)

    if not chr_names:
        raise ValueError(
            f"No main chromosomes matching ^chr[0-9]+$ were found in {chr_bed}"
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


def filter_main_chr_rows(
    rows: list[dict[str, str]],
    valid_chr: set[str],
) -> list[dict[str, str]]:
    kept = []

    for row in rows:
        chr1 = (row.get("chr1") or "").strip()
        chr2 = (row.get("chr2") or "").strip()

        if chr1 in valid_chr and chr2 in valid_chr:
            kept.append(row)

    return kept


def write_links_file(rows: list[dict[str, str]], out_path: Path) -> None:
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
    species: str,
    out_conf: Path,
    karyotype_path: Path,
    full_links_path: Path,
) -> None:
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
label_font       = helvetica
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
file          = {full_links_path.name}
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
file              = circos_{species}
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


def main() -> None:
    args = parse_args()

    species = args.species
    circos_tsv = Path(args.circos_links)
    chr_bed = Path(args.chr_bed)
    species_dir = Path(args.output_dir)

    links_by_branch_dir = species_dir / "links_by_redip_branch"
    species_dir.mkdir(parents=True, exist_ok=True)
    links_by_branch_dir.mkdir(parents=True, exist_ok=True)

    if not chr_bed.exists():
        raise FileNotFoundError(f"Missing chromosome BED for {species}: {chr_bed}")

    if not circos_tsv.exists():
        raise FileNotFoundError(f"Missing Circos links TSV for {species}: {circos_tsv}")

    karyotype_path = species_dir / "karyotype.txt"
    full_links_path = species_dir / f"{species}.links.tsv"
    conf_path = species_dir / "circos.conf"

    chr_names = make_karyotype(chr_bed, karyotype_path)
    valid_chr = set(chr_names)

    all_rows = load_links(circos_tsv)
    kept_rows = filter_main_chr_rows(all_rows, valid_chr)

    write_links_file(kept_rows, full_links_path)

    branch_rows = defaultdict(list)
    for row in kept_rows:
        branch = (row.get("redip_branch") or "").strip()
        if branch:
            branch_rows[branch].append(row)

    for branch in sorted(branch_rows, key=branch_sort_key):
        out_branch = links_by_branch_dir / f"{species}.branch_{branch}.links.tsv"
        write_links_file(branch_rows[branch], out_branch)

    write_circos_conf(
        species=species,
        out_conf=conf_path,
        karyotype_path=karyotype_path,
        full_links_path=full_links_path,
    )

    summary_path = species_dir / "README.txt"
    with open(summary_path, "w") as out:
        out.write(f"Species: {species}\n")
        out.write(f"Karyotype: {karyotype_path.name}\n")
        out.write(f"Full links: {full_links_path.name}\n")
        out.write(f"Circos config: {conf_path.name}\n")
        out.write(f"Rows retained on main chromosomes: {len(kept_rows)}\n")
        out.write("Per-branch files:\n")

        for branch in sorted(branch_rows, key=branch_sort_key):
            out.write(
                f"  links_by_redip_branch/{species}.branch_{branch}.links.tsv"
                f"\t{len(branch_rows[branch])} rows\n"
            )

    print(f"Prepared Circos inputs for {species}")
    print(f"Output directory: {species_dir}")


if __name__ == "__main__":
    main()
