#!/usr/bin/env python3

# Adapted from Emily Haley's 2026-04 filterTandemDups.py

import pandas as pd
import argparse
from pathlib import Path


def parse_bool(value: str) -> bool:
    v = str(value).strip().lower()
    if v in {"yes", "true", "1"}:
        return True
    if v in {"no", "false", "0", ""}:
        return False
    raise SystemExit(
        f"Could not parse boolean value '{value}'. Use yes/no, true/false, or 1/0."
    )


def load_outgroup_genomes(genomes_tsv: Path, require_outgroup: bool) -> set[str]:
    df = pd.read_csv(genomes_tsv, sep="\t", dtype=str).fillna("")

    if "genome_id" not in df.columns:
        raise SystemExit("genomes TSV must contain a 'genome_id' column")

    if "outgroup" not in df.columns:
        if require_outgroup:
            raise SystemExit(
                "The genomes TSV must contain an 'outgroup' column when --require-outgroup is used."
            )
        return set()

    df["genome_id"] = df["genome_id"].astype(str).str.strip()
    df["outgroup"] = df["outgroup"].astype(str).str.strip()

    return {
        row["genome_id"]
        for _, row in df.iterrows()
        if parse_bool(row["outgroup"])
    }


def main():
    parser = argparse.ArgumentParser(
        description="Collapse tandem duplicates with configurable ord gap"
    )
    parser.add_argument("--infile", required=True, help="Input PASS pangenes TSV")
    parser.add_argument("--genomes-tsv", required=True, help="Input genomes TSV")
    parser.add_argument("--outfile_filtered", required=True, help="Filtered output TSV")
    parser.add_argument("--outfile_tandems", required=True, help="Tandem report TSV")
    parser.add_argument("--outfile_og_list", required=True, help="Rewritten OG list")
    parser.add_argument(
        "--max_ord_gap",
        type=int,
        default=1,
        help="Maximum allowed gap in gene order to define tandem cluster (default: 1)"
    )
    parser.add_argument(
        "--require-outgroup",
        action="store_true",
        help="Require each retained OG to include at least one outgroup genome"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.infile, sep="\t", dtype=str)

    required_cols = ["og", "genome", "chr", "id", "ord"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in input TSV: {', '.join(missing)}")

    for col in required_cols:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df["ord_num"] = pd.to_numeric(df["ord"], errors="coerce")

    kept_rows = []
    tandem_records = []

    for (og, genome, chr_id), subdf in df.groupby(["og", "genome", "chr"], sort=False):
        subdf = subdf.sort_values("ord_num")

        if len(subdf) == 1:
            kept_rows.append(subdf.iloc[0])
            continue

        cluster = [subdf.iloc[0]]

        for i in range(1, len(subdf)):
            prev = cluster[-1]
            curr = subdf.iloc[i]

            if (
                pd.notna(prev["ord_num"]) and
                pd.notna(curr["ord_num"]) and
                curr["ord_num"] <= prev["ord_num"] + args.max_ord_gap
            ):
                cluster.append(curr)
            else:
                kept = cluster[0]
                kept_rows.append(kept)

                if len(cluster) > 1:
                    removed = [row["id"] for row in cluster[1:]]
                    tandem_records.append({
                        "og": og,
                        "genome": genome,
                        "chr": chr_id,
                        "kept_gene": kept["id"],
                        "removed_genes": ",".join(removed)
                    })

                cluster = [curr]

        kept = cluster[0]
        kept_rows.append(kept)

        if len(cluster) > 1:
            removed = [row["id"] for row in cluster[1:]]
            tandem_records.append({
                "og": og,
                "genome": genome,
                "chr": chr_id,
                "kept_gene": kept["id"],
                "removed_genes": ",".join(removed)
            })

    out_df = pd.DataFrame(kept_rows).drop(columns=["ord_num"], errors="ignore")
    tandem_df = pd.DataFrame(tandem_records)

    outgroup_genomes = load_outgroup_genomes(Path(args.genomes_tsv), args.require_outgroup)

    og_stats = out_df.groupby("og").agg(
        n_rows=("id", "size"),
        n_genomes=("genome", "nunique")
    ).reset_index()

    if args.require_outgroup:
        outgroup_flag = (
            out_df.assign(is_outgroup=out_df["genome"].isin(outgroup_genomes))
            .groupby("og")["is_outgroup"]
            .any()
            .reset_index()
        )
        og_stats = og_stats.merge(outgroup_flag, on="og", how="left")
        og_stats["is_outgroup"] = og_stats["is_outgroup"].fillna(False)

        keep_ogs = set(
            og_stats.loc[
                (og_stats["n_rows"] >= 4) &
                (og_stats["n_genomes"] >= 4) &
                (og_stats["is_outgroup"]),
                "og"
            ]
        )
    else:
        keep_ogs = set(
            og_stats.loc[
                (og_stats["n_rows"] >= 4) &
                (og_stats["n_genomes"] >= 4),
                "og"
            ]
        )

    out_df = out_df[out_df["og"].isin(keep_ogs)].copy()

    out_df.to_csv(args.outfile_filtered, sep="\t", index=False)

    with open(args.outfile_og_list, "w") as f:
        for og in sorted(keep_ogs):
            f.write(f"{og}\n")

    if len(tandem_df) == 0:
        print("No tandem clusters found.")
        pd.DataFrame(columns=["og", "genome", "chr", "kept_gene", "removed_genes"]) \
            .to_csv(args.outfile_tandems, sep="\t", index=False)
    else:
        tandem_df.to_csv(args.outfile_tandems, sep="\t", index=False)

    print("--------------------------------------------------")
    print(f"Input rows: {len(df)}")
    print(f"Output rows: {len(out_df)}")
    print(f"Removed tandem genes: {len(df) - len(out_df)}")
    print(f"Tandem clusters found: {len(tandem_df)}")
    print(f"Filtered file: {args.outfile_filtered}")
    print(f"Tandem report: {args.outfile_tandems}")
    print(f"Rewritten OG list: {args.outfile_og_list}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    main()
