#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


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
        description="Collapse tandem duplicates, rewrite the filtered TSV, and regenerate the OG list."
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
        help="Maximum allowed difference in ord values to define a tandem cluster (default: 1)"
    )
    parser.add_argument(
        "--require-outgroup",
        action="store_true",
        help="Require each retained OG to include at least one outgroup genome"
    )
    args = parser.parse_args()

    infile = Path(args.infile)
    genomes_tsv = Path(args.genomes_tsv)
    outfile_filtered = Path(args.outfile_filtered)
    outfile_tandems = Path(args.outfile_tandems)
    outfile_og_list = Path(args.outfile_og_list)

    df = pd.read_csv(infile, sep="\t", dtype=str)

    required_cols = ["og", "genome", "chr", "id", "ord"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in input TSV: {', '.join(missing)}")

    for col in required_cols:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df["ord_num"] = pd.to_numeric(df["ord"], errors="coerce")

    input_rows = len(df)
    input_ogs = df["og"].nunique()

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

            same_cluster = (
                pd.notna(prev["ord_num"]) and
                pd.notna(curr["ord_num"]) and
                curr["ord_num"] <= prev["ord_num"] + args.max_ord_gap
            )

            if same_cluster:
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
                        "removed_genes": ",".join(removed),
                        "cluster_size": len(cluster),
                        "n_removed": len(cluster) - 1
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
                "removed_genes": ",".join(removed),
                "cluster_size": len(cluster),
                "n_removed": len(cluster) - 1
            })

    collapsed_df = pd.DataFrame(kept_rows).drop(columns=["ord_num"], errors="ignore")
    tandem_df = pd.DataFrame(tandem_records)

    rows_after_collapse = len(collapsed_df)
    ogs_after_collapse = collapsed_df["og"].nunique() if not collapsed_df.empty else 0
    tandem_clusters_found = len(tandem_df)
    rows_removed_by_tandem = input_rows - rows_after_collapse

    outgroup_genomes = load_outgroup_genomes(genomes_tsv, args.require_outgroup)

    if collapsed_df.empty:
        raise SystemExit("No rows remain after tandem collapsing.")

    og_stats = collapsed_df.groupby("og").agg(
        n_rows=("id", "size"),
        n_genomes=("genome", "nunique")
    ).reset_index()

    if args.require_outgroup:
        outgroup_flag = (
            collapsed_df.assign(is_outgroup=collapsed_df["genome"].isin(outgroup_genomes))
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

    final_df = collapsed_df[collapsed_df["og"].isin(keep_ogs)].copy()

    final_rows = len(final_df)
    final_ogs = final_df["og"].nunique() if not final_df.empty else 0
    rows_removed_by_og_refilter = rows_after_collapse - final_rows

    outfile_filtered.parent.mkdir(parents=True, exist_ok=True)
    outfile_tandems.parent.mkdir(parents=True, exist_ok=True)
    outfile_og_list.parent.mkdir(parents=True, exist_ok=True)

    if final_df.empty:
        raise SystemExit(
            "No rows remain after re-filtering OGs following tandem collapse. "
            "No output files were written."
        )

    final_df.to_csv(outfile_filtered, sep="\t", index=False)

    with open(outfile_og_list, "w") as f:
        for og in sorted(keep_ogs):
            f.write(f"{og}\n")

    if tandem_df.empty:
        pd.DataFrame(
            columns=["og", "genome", "chr", "kept_gene", "removed_genes", "cluster_size", "n_removed"]
        ).to_csv(outfile_tandems, sep="\t", index=False)
    else:
        tandem_df.to_csv(outfile_tandems, sep="\t", index=False)

    print("--------------------------------------------------")
    print("collapse_tandems summary")
    print("--------------------------------------------------")
    print(f"Input file: {infile}")
    print(f"Input rows: {input_rows}")
    print(f"Input OGs: {input_ogs}")
    print(f"Max ord gap used: {args.max_ord_gap}")
    print(f"Outgroup requirement enabled: {'yes' if args.require_outgroup else 'no'}")
    if args.require_outgroup:
        print(f"Outgroup genomes loaded: {len(outgroup_genomes)}")
    print("--------------------------------------------------")
    print(f"Rows after tandem collapse: {rows_after_collapse}")
    print(f"OGs after tandem collapse: {ogs_after_collapse}")
    print(f"Tandem clusters found: {tandem_clusters_found}")
    print(f"Rows removed by tandem collapse: {rows_removed_by_tandem}")
    print("--------------------------------------------------")
    print(f"Final rows after OG re-filtering: {final_rows}")
    print(f"Final OGs retained: {final_ogs}")
    print(f"Rows dropped by OG re-filtering after collapse: {rows_removed_by_og_refilter}")
    print("--------------------------------------------------")
    print(f"Filtered TSV written to: {outfile_filtered}")
    print(f"Tandem report written to: {outfile_tandems}")
    print(f"OG list written to: {outfile_og_list}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    main()
