#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd

TRAILING_DOTNUM = re.compile(r"\.\d+$")

def norm_drop_dotnum(s: str) -> str:
    s = s.strip()
    return TRAILING_DOTNUM.sub("", s)

def fasta_iter(path: Path):
    header = None
    seq_chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)

def build_cds_index_for_genome(cds_path: Path, source: str):
    """
    Returns a dict mapping multiple possible keys -> sequence,
    plus the raw headers list for fallback checks.
    """
    index = {}
    headers_first_field = []

    for hdr, seq in fasta_iter(cds_path):
        tokens = hdr.split()
        first = tokens[0] if tokens else hdr
        headers_first_field.append(first)

        # Always index exact first field and full header as basic anchors
        index[first] = seq
        index[hdr] = seq

        if source == "ensembl":
            # Field 4 (space-separated), strip "gene" text, drop trailing .#
            if len(tokens) >= 4:
                t = tokens[3]
                t = t.replace("gene", "").replace(":", "").strip()
                t = norm_drop_dotnum(t)
                if t:
                    index[t] = seq

        elif source == "phytozome":
            # Index first field normalized without trailing .#
            index[norm_drop_dotnum(first)] = seq

            # Also index field 5 with ID= stripped, normalized
            if len(tokens) >= 5:
                t = tokens[4]
                if t.startswith("ID="):
                    t = t[3:]
                t = norm_drop_dotnum(t)
                if t:
                    index[t] = seq

        elif source == "ncbi":
            # NCBI: we’ll try exact match first, but also index normalized first field
            index[norm_drop_dotnum(first)] = seq

    return index, headers_first_field

def find_seq_for_id(genome_source: str, cds_index: dict, headers_first_field: list, target_id: str):
    target_id = str(target_id).strip()

    if genome_source == "ensembl":
        # We expect we already indexed the extracted gene token (field 4 processed),
        # but also allow direct lookup by id
        if target_id in cds_index:
            return cds_index[target_id]
        t = norm_drop_dotnum(target_id)
        if t in cds_index:
            return cds_index[t]
        return None

    if genome_source == "ncbi":
        # Step 1: exact match between id and cds fasta header
        # Interpret as: exact match against first field OR full header key
        if target_id in cds_index:
            return cds_index[target_id]

        # Step 2: match somewhere within the first header field
        for first in headers_first_field:
            if target_id in first:
                seq = cds_index.get(first)
                if seq:
                    return seq
        return None

    if genome_source == "phytozome":
        # Step 1: match between id and header ignoring trailing .# in header
        # We indexed norm_drop_dotnum(first_field)
        if target_id in cds_index:
            return cds_index[target_id]
        if norm_drop_dotnum(target_id) in cds_index:
            return cds_index[norm_drop_dotnum(target_id)]

        # Step 2: field 5 ID=... normalized was indexed too
        # so if step 1 failed, we still try direct lookup
        # (already covered above), plus try scanning for ID=
        for key in (target_id, norm_drop_dotnum(target_id)):
            if key in cds_index:
                return cds_index[key]
        return None

    # Default: direct
    return cds_index.get(target_id)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenes-pass", required=True)
    ap.add_argument("--genomes-tsv", required=True)
    ap.add_argument("--cds-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--og-list", required=True)
    args = ap.parse_args()

    ppass = pd.read_csv(args.pangenes_pass, sep="\t", dtype=str)
    genomes = pd.read_csv(args.genomes_tsv, sep="\t", dtype=str).rename(
        columns={"genome_id": "genome", "genome_source": "source"}
    )
    genomes["genome"] = genomes["genome"].str.strip()
    genomes["source"] = genomes["source"].str.lower().str.strip()
    src = dict(zip(genomes["genome"], genomes["source"]))

    cols = list(ppass.columns)
    def col_by_name_or_pos(name: str, pos_1based: int) -> str:
        if name in ppass.columns:
            return name
        idx = pos_1based - 1
        if idx < 0 or idx >= len(cols):
            raise SystemExit(f"Expected column {pos_1based} for '{name}' but file has only {len(cols)} columns")
        return cols[idx]

    col_og     = col_by_name_or_pos("og", 7)
    col_genome = col_by_name_or_pos("genome", 6)
    col_id     = col_by_name_or_pos("id", 9)
    col_chr    = col_by_name_or_pos("chr", 10)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Build CDS indices per genome only once
    cds_dir = Path(args.cds_dir)
    cds_cache = {}

    def get_cds_index(genome_id: str):
        if genome_id in cds_cache:
            return cds_cache[genome_id]
        cds_path = cds_dir / f"{genome_id}.cds"
        if not cds_path.exists():
            raise SystemExit(f"Missing CDS fasta: {cds_path}")
        genome_source = src.get(genome_id)
        if genome_source is None:
            raise SystemExit(f"Genome {genome_id} not found in genomes.tsv")
        idx, first_fields = build_cds_index_for_genome(cds_path, genome_source)
        cds_cache[genome_id] = (genome_source, idx, first_fields)
        return cds_cache[genome_id]

    ogs = sorted(ppass[col_og].astype(str).unique().tolist())
    with open(args.og_list, "w") as f:
        for og in ogs:
            f.write(f"{og}\n")

    # Write one fasta per OG
    for og in ogs:
        sub = ppass[ppass[col_og].astype(str) == og].copy()

        out_fa = outdir / f"og_{og}.fasta"
        wrote_any = False

        with open(out_fa, "w") as out:
            for _, row in sub.iterrows():
                genome_id = str(row[col_genome]).strip()
                rec_id = str(row[col_id]).strip()
                chr_id = str(row[col_chr]).strip()

                genome_source, cds_index, first_fields = get_cds_index(genome_id)
                seq = find_seq_for_id(genome_source, cds_index, first_fields, rec_id)
                if seq is None:
                    raise SystemExit(
                        f"Failed to find CDS for genome={genome_id} source={genome_source} id={rec_id} og={og}"
                    )

                out_hdr = f"{genome_id}|{chr_id}|{rec_id}"
                out.write(f">{out_hdr}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
                wrote_any = True

        if not wrote_any:
            raise SystemExit(f"No sequences written for OG {og} (unexpected)")

    # Completion marker is created by Snakemake rule (og_fastas.done)

if __name__ == "__main__":
    main()
