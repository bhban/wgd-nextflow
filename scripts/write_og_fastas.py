#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


TRAILING_DOTNUM = re.compile(r"\.\d+$")


def norm_drop_dotnum(s: str) -> str:
    s = str(s).strip()
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


def add_index_key(index: dict, key: str, seq: str):
    key = str(key).strip()
    if key:
        index[key] = seq


def build_sequence_index_for_genome(fasta_path: Path, source: str):
    """
    Returns:
      index: dict mapping possible gene/transcript/header keys to sequence
      headers_first_field: list of first FASTA header fields for fallback matching

    This is used for both CDS and protein FASTA files. The source-specific
    matching rules are intentionally inherited from the original CDS-only script.
    """
    index = {}
    headers_first_field = []

    for hdr, seq in fasta_iter(fasta_path):
        tokens = hdr.split()
        first = tokens[0] if tokens else hdr
        headers_first_field.append(first)

        # Basic anchors
        add_index_key(index, first, seq)
        add_index_key(index, hdr, seq)
        add_index_key(index, norm_drop_dotnum(first), seq)

        if source == "ensembl":
            # Original logic:
            # Field 4, strip "gene" text and ":", then drop trailing .#
            if len(tokens) >= 4:
                t = tokens[3]
                t = t.replace("gene", "").replace(":", "").strip()
                t = norm_drop_dotnum(t)
                add_index_key(index, t, seq)

            # Extra conservative anchors, useful if headers differ slightly
            for t in tokens:
                clean = t.replace("gene", "").replace(":", "").strip()
                add_index_key(index, clean, seq)
                add_index_key(index, norm_drop_dotnum(clean), seq)

        elif source == "phytozome":
            # Original logic:
            # Index first field normalized without trailing .#
            add_index_key(index, norm_drop_dotnum(first), seq)

            # Original logic:
            # Also index field 5 with ID= stripped, normalized
            if len(tokens) >= 5:
                t = tokens[4]
                if t.startswith("ID="):
                    t = t[3:]
                t = norm_drop_dotnum(t)
                add_index_key(index, t, seq)

            # Extra conservative ID= scan across tokens
            for t in tokens:
                if t.startswith("ID="):
                    t = t[3:]
                    add_index_key(index, t, seq)
                    add_index_key(index, norm_drop_dotnum(t), seq)

        elif source == "ncbi":
            # Original logic:
            # Try exact match first, but also index normalized first field
            add_index_key(index, norm_drop_dotnum(first), seq)

            # Extra conservative anchors for common NCBI-ish key-value tokens
            for t in tokens:
                clean = t.strip()
                add_index_key(index, clean, seq)
                add_index_key(index, norm_drop_dotnum(clean), seq)

                for prefix in ("gene=", "GeneID:", "protein_id=", "transcript_id=", "locus_tag="):
                    if clean.startswith(prefix):
                        val = clean.replace(prefix, "", 1).strip()
                        add_index_key(index, val, seq)
                        add_index_key(index, norm_drop_dotnum(val), seq)

    return index, headers_first_field


def find_seq_for_id(genome_source: str, seq_index: dict, headers_first_field: list, target_id: str):
    target_id = str(target_id).strip()
    target_norm = norm_drop_dotnum(target_id)

    if target_id in seq_index:
        return seq_index[target_id]

    if target_norm in seq_index:
        return seq_index[target_norm]

    if genome_source == "ensembl":
        return None

    if genome_source == "ncbi":
        # Original fallback:
        # match target ID somewhere within the first FASTA header field
        for first in headers_first_field:
            if target_id in first or target_norm in norm_drop_dotnum(first):
                seq = seq_index.get(first) or seq_index.get(norm_drop_dotnum(first))
                if seq:
                    return seq
        return None

    if genome_source == "phytozome":
        return None

    return None


def find_fasta_for_genome(seq_dir: Path, genome_id: str, kind: str):
    """
    Finds a genome FASTA in the staged working directory.

    The preferred filenames are:
      <genome_id>.cds for CDS
      <genome_id>.pep for protein

    The fallback accepts common FASTA suffixes because different parts of the
    pipeline may preserve original extensions.
    """
    if kind == "cds":
        preferred_names = [
            f"{genome_id}.cds",
            f"{genome_id}.cds.fa",
            f"{genome_id}.cds.fasta",
            f"{genome_id}.fna",
            f"{genome_id}.fa",
            f"{genome_id}.fasta",
        ]
    elif kind == "protein":
        preferred_names = [
            f"{genome_id}.pep",
            f"{genome_id}.pep.fa",
            f"{genome_id}.pep.fasta",
            f"{genome_id}.faa",
            f"{genome_id}.fa",
            f"{genome_id}.fasta",
        ]
    else:
        raise ValueError(f"Unknown sequence kind: {kind}")

    for name in preferred_names:
        p = seq_dir / name
        if p.exists():
            return p

    matches = sorted(seq_dir.glob(f"{genome_id}.*"))
    matches = [p for p in matches if p.is_file()]

    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:
        names = ", ".join(str(p.name) for p in matches)
        raise SystemExit(
            f"Multiple possible {kind} FASTA files found for genome={genome_id} in {seq_dir}: {names}"
        )

    raise SystemExit(f"Missing {kind} FASTA for genome={genome_id} in {seq_dir}")


def load_genomes(genomes_tsv: str):
    genomes = pd.read_csv(genomes_tsv, sep="\t", dtype=str).rename(
        columns={
            "genome_id": "genome",
            "genome_source": "source",
        }
    )

    if "genome" not in genomes.columns:
        raise SystemExit("genomes.tsv must contain a genome_id column, or a genome column")

    if "source" not in genomes.columns:
        raise SystemExit("genomes.tsv must contain a genome_source column, or a source column")

    genomes["genome"] = genomes["genome"].astype(str).str.strip()
    genomes["source"] = genomes["source"].astype(str).str.lower().str.strip()

    return dict(zip(genomes["genome"], genomes["source"]))


def col_by_name_or_pos(df: pd.DataFrame, name: str, pos_1based: int) -> str:
    if name in df.columns:
        return name

    cols = list(df.columns)
    idx = pos_1based - 1

    if idx < 0 or idx >= len(cols):
        raise SystemExit(
            f"Expected column {pos_1based} for '{name}' but file has only {len(cols)} columns"
        )

    return cols[idx]


def sort_ogs(ogs):
    ogs = [str(x).strip() for x in ogs]

    def key(x):
        try:
            return (0, int(x))
        except ValueError:
            return (1, x)

    return sorted(ogs, key=key)


def write_wrapped_fasta_record(out, header: str, seq: str, width: int = 60):
    out.write(f">{header}\n")
    for i in range(0, len(seq), width):
        out.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenes-pass", required=True)
    ap.add_argument("--genomes-tsv", required=True)

    ap.add_argument("--cds-dir", required=True)
    ap.add_argument("--protein-dir", required=True)

    ap.add_argument("--outdir-cds", required=True)
    ap.add_argument("--outdir-aa", required=True)

    ap.add_argument("--og-list", required=True)

    args = ap.parse_args()

    ppass = pd.read_csv(args.pangenes_pass, sep="\t", dtype=str)
    src = load_genomes(args.genomes_tsv)

    col_og = col_by_name_or_pos(ppass, "og", 7)
    col_genome = col_by_name_or_pos(ppass, "genome", 6)
    col_id = col_by_name_or_pos(ppass, "id", 9)
    col_chr = col_by_name_or_pos(ppass, "chr", 10)

    cds_dir = Path(args.cds_dir)
    protein_dir = Path(args.protein_dir)

    outdir_cds = Path(args.outdir_cds)
    outdir_aa = Path(args.outdir_aa)

    outdir_cds.mkdir(parents=True, exist_ok=True)
    outdir_aa.mkdir(parents=True, exist_ok=True)

    cds_cache = {}
    protein_cache = {}

    def get_sequence_index(genome_id: str, kind: str):
        genome_id = str(genome_id).strip()

        if genome_id not in src:
            raise SystemExit(f"Genome {genome_id} not found in genomes.tsv")

        genome_source = src[genome_id]

        if kind == "cds":
            cache = cds_cache
            seq_dir = cds_dir
        elif kind == "protein":
            cache = protein_cache
            seq_dir = protein_dir
        else:
            raise ValueError(f"Unknown sequence kind: {kind}")

        if genome_id in cache:
            return cache[genome_id]

        seq_path = find_fasta_for_genome(seq_dir, genome_id, kind)
        idx, first_fields = build_sequence_index_for_genome(seq_path, genome_source)

        cache[genome_id] = {
            "source": genome_source,
            "path": seq_path,
            "index": idx,
            "first_fields": first_fields,
        }

        return cache[genome_id]

    ogs = sort_ogs(ppass[col_og].astype(str).unique().tolist())

    with open(args.og_list, "w") as f:
        for og in ogs:
            f.write(f"{og}\n")

    for og in ogs:
        sub = ppass[ppass[col_og].astype(str) == og].copy()

        out_cds_fa = outdir_cds / f"og_{og}.fasta"
        out_aa_fa = outdir_aa / f"og_{og}.fasta"

        wrote_cds = False
        wrote_aa = False

        with open(out_cds_fa, "w") as out_cds, open(out_aa_fa, "w") as out_aa:
            for _, row in sub.iterrows():
                genome_id = str(row[col_genome]).strip()
                rec_id = str(row[col_id]).strip()
                chr_id = str(row[col_chr]).strip()

                out_hdr = f"{genome_id}|{chr_id}|{rec_id}"

                cds_info = get_sequence_index(genome_id, "cds")
                cds_seq = find_seq_for_id(
                    cds_info["source"],
                    cds_info["index"],
                    cds_info["first_fields"],
                    rec_id,
                )

                if cds_seq is None:
                    raise SystemExit(
                        "Failed to find CDS sequence "
                        f"for genome={genome_id} "
                        f"source={cds_info['source']} "
                        f"id={rec_id} "
                        f"og={og} "
                        f"file={cds_info['path']}"
                    )

                protein_info = get_sequence_index(genome_id, "protein")
                aa_seq = find_seq_for_id(
                    protein_info["source"],
                    protein_info["index"],
                    protein_info["first_fields"],
                    rec_id,
                )

                if aa_seq is None:
                    raise SystemExit(
                        "Failed to find protein sequence "
                        f"for genome={genome_id} "
                        f"source={protein_info['source']} "
                        f"id={rec_id} "
                        f"og={og} "
                        f"file={protein_info['path']}"
                    )

                write_wrapped_fasta_record(out_cds, out_hdr, cds_seq)
                wrote_cds = True

                write_wrapped_fasta_record(out_aa, out_hdr, aa_seq)
                wrote_aa = True

        if not wrote_cds:
            raise SystemExit(f"No CDS sequences written for OG {og}")

        if not wrote_aa:
            raise SystemExit(f"No protein sequences written for OG {og}")


if __name__ == "__main__":
    main()
