#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
from collections import defaultdict
import pandas as pd


def open_text_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def resolve_columns(header_cols):
    cols = list(header_cols)

    def col_by_name_or_pos(name: str, pos_1based: int) -> str:
        if name in cols:
            return name
        idx = pos_1based - 1
        if idx < 0 or idx >= len(cols):
            raise SystemExit(
                f"Expected column {pos_1based} for '{name}' but file has only {len(cols)} columns"
            )
        return cols[idx]

    col_flag = col_by_name_or_pos("flag", 8)
    col_id = col_by_name_or_pos("id", 9)
    col_og = col_by_name_or_pos("og", 7)
    col_genome = col_by_name_or_pos("genome", 6)

    return col_flag, col_id, col_og, col_genome


def first_chunk_and_columns(path: Path, chunksize: int):
    with open_text_maybe_gzip(path) as f:
        reader = pd.read_csv(f, sep="\t", dtype=str, chunksize=chunksize)
        first_chunk = next(reader, None)
    if first_chunk is None:
        raise SystemExit(f"File is empty: {path}")
    return first_chunk, list(first_chunk.columns)


def iter_chunks(path: Path, chunksize: int):
    with open_text_maybe_gzip(path) as f:
        yield from pd.read_csv(f, sep="\t", dtype=str, chunksize=chunksize)


def normalize_key_columns(df: pd.DataFrame, col_flag: str, col_id: str, col_og: str, col_genome: str) -> pd.DataFrame:
    # Only normalize the columns needed for filtering/grouping.
    # Fill NA before string ops so behavior is stable.
    df[col_flag] = df[col_flag].fillna("").astype(str).str.strip()
    df[col_id] = df[col_id].fillna("").astype(str).str.strip()
    df[col_og] = df[col_og].fillna("").astype(str).str.strip()
    df[col_genome] = df[col_genome].fillna("").astype(str).str.strip()
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genespace-wd", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-og-list", required=True)
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()

    wd = Path(args.genespace_wd)
    pang_dir = wd / "pangenes"
    if not pang_dir.exists():
        raise SystemExit(f"Missing pangenes dir: {pang_dir}")

    files = sorted(pang_dir.glob("*_pangenes.txt.gz"))
    if not files:
        raise SystemExit(f"No pangenes files found in: {pang_dir}")

    # Determine columns from the first file header/chunk
    first_chunk, first_cols = first_chunk_and_columns(files[0], args.chunksize)
    col_flag, col_id, col_og, col_genome = resolve_columns(first_cols)

    # Pass 1: determine which OGs pass after global PASS filter + global first-id dedup
    seen_ids = set()
    row_count = defaultdict(int)
    genomes_by_og = defaultdict(set)

    for fp in files:
        for chunk in iter_chunks(fp, args.chunksize):
            # Ensure expected columns exist
            missing = [c for c in [col_flag, col_id, col_og, col_genome] if c not in chunk.columns]
            if missing:
                raise SystemExit(f"Missing expected columns {missing} in file: {fp}")

            chunk = normalize_key_columns(chunk, col_flag, col_id, col_og, col_genome)

            # PASS only
            chunk = chunk[chunk[col_flag] == "PASS"]
            if chunk.empty:
                continue

            # Preserve original behavior: keep first occurrence of each id in global order
            keep_mask = []
            for gene_id in chunk[col_id]:
                if gene_id in seen_ids:
                    keep_mask.append(False)
                else:
                    seen_ids.add(gene_id)
                    keep_mask.append(True)

            chunk = chunk.loc[keep_mask]
            if chunk.empty:
                continue

            # Update OG stats
            for og, subdf in chunk.groupby(col_og, sort=False, dropna=False):
                row_count[og] += len(subdf)
                genomes_by_og[og].update(subdf[col_genome].tolist())

    keep_ogs = {
        og for og in row_count
        if row_count[og] >= 4 and len(genomes_by_og[og]) >= 4
    }

    # Write OG list
    out_og = Path(args.out_og_list)
    out_og.parent.mkdir(parents=True, exist_ok=True)
    with open(out_og, "w") as f:
        for og in sorted(keep_ogs):
            f.write(f"{og}\n")

    # Pass 2: write all retained rows with all original columns
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    seen_ids = set()
    wrote_header = False

    for fp in files:
        for chunk in iter_chunks(fp, args.chunksize):
            missing = [c for c in [col_flag, col_id, col_og, col_genome] if c not in chunk.columns]
            if missing:
                raise SystemExit(f"Missing expected columns {missing} in file: {fp}")

            chunk = normalize_key_columns(chunk, col_flag, col_id, col_og, col_genome)

            # PASS only
            chunk = chunk[chunk[col_flag] == "PASS"]
            if chunk.empty:
                continue

            # Global first-occurrence dedup by id, same as original script
            keep_mask = []
            for gene_id in chunk[col_id]:
                if gene_id in seen_ids:
                    keep_mask.append(False)
                else:
                    seen_ids.add(gene_id)
                    keep_mask.append(True)

            chunk = chunk.loc[keep_mask]
            if chunk.empty:
                continue

            # Keep only OGs that pass thresholds
            chunk = chunk[chunk[col_og].isin(keep_ogs)]
            if chunk.empty:
                continue

            chunk.to_csv(
                out_tsv,
                sep="\t",
                index=False,
                mode="a",
                header=not wrote_header
            )
            wrote_header = True

    if not wrote_header:
        raise SystemExit("No rows passed filtering; output TSV was not written.")

    if out_tsv.stat().st_size == 0:
        raise SystemExit(f"Empty output TSV: {out_tsv}")
    if out_og.stat().st_size == 0:
        raise SystemExit(f"Empty OG list: {out_og}")


if __name__ == "__main__":
    main()
