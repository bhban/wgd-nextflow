#!/usr/bin/env python3

# Emily Haley's 2026-04 filterTandemDups.py

import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Collapse tandem duplicates with configurable ord gap"
    )
    parser.add_argument("--infile", required=True, help="Input PASS pangenes TSV")
    parser.add_argument("--outfile_filtered", required=True, help="Filtered output TSV")
    parser.add_argument("--outfile_tandems", required=True, help="Tandem report TSV")
    parser.add_argument("--max_ord_gap", type=int, default=1,
                        help="Maximum allowed gap in gene order to define tandem cluster (default: 1)")
    args = parser.parse_args()

    df = pd.read_csv(args.infile, sep="\t", dtype=str)

    # Clean columns
    for col in ["og", "genome", "chr", "id", "ord"]:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df["ord_num"] = pd.to_numeric(df["ord"], errors="coerce")

    kept_rows = []
    tandem_records = []

    # Group by OG, genome, chromosome
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
                # Same tandem cluster
                cluster.append(curr)
            else:
                # Resolve cluster
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

        # Final cluster
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

    # Create outputs
    out_df = pd.DataFrame(kept_rows).drop(columns=["ord_num"], errors="ignore")
    tandem_df = pd.DataFrame(tandem_records)

    # Write files
    out_df.to_csv(args.outfile_filtered, sep="\t", index=False)

    if len(tandem_df) == 0:
        print("No tandem clusters found.")
        pd.DataFrame(columns=["og", "genome", "chr", "kept_gene", "removed_genes"])\
            .to_csv(args.outfile_tandems, sep="\t", index=False)
    else:
        tandem_df.to_csv(args.outfile_tandems, sep="\t", index=False)

    # Summary
    print("--------------------------------------------------")
    print(f"Input rows: {len(df)}")
    print(f"Output rows: {len(out_df)}")
    print(f"Removed tandem genes: {len(df) - len(out_df)}")
    print(f"Tandem clusters found: {len(tandem_df)}")
    print(f"Filtered file: {args.outfile_filtered}")
    print(f"Tandem report: {args.outfile_tandems}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    main()
