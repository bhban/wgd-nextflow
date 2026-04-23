#!/usr/bin/env python3

# Adapted from Emily Haley's filterTandemDups.py
# Now uses  arrayID and isArrayRepfrom GENESPACE's results/combBed.txt to collapse tandem arrays if possible

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


def normalize_bool_str(value: str) -> str:
    v = str(value).strip().upper()
    if v in {"TRUE", "T", "1", "YES"}:
        return "TRUE"
    if v in {"FALSE", "F", "0", "NO"}:
        return "FALSE"
    return ""


def choose_kept_row(cluster: pd.DataFrame) -> tuple[pd.Series, str, int]:
    """
    Choose the representative row for one collapsed tandem cluster.

    Priority:
    1. Keep a row with isArrayRep == TRUE
    2. If multiple TRUE rows, choose the one with lowest ord_num then id
    3. If no TRUE rows, fall back to the row with lowest ord_num then id

    Returns:
        kept_row, selection_method, n_true_reps
    """
    cluster = cluster.sort_values(["ord_num", "id"], na_position="last").copy()
    true_rep_mask = cluster["isArrayRep_norm"] == "TRUE"
    true_reps = cluster.loc[true_rep_mask].copy()
    n_true_reps = len(true_reps)

    if n_true_reps >= 1:
        kept = true_reps.sort_values(["ord_num", "id"], na_position="last").iloc[0]
        if n_true_reps == 1:
            selection_method = "isArrayRep_TRUE"
        else:
            selection_method = "multiple_TRUE_lowest_ord"
    else:
        kept = cluster.iloc[0]
        selection_method = "fallback_lowest_ord"

    return kept, selection_method, n_true_reps


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Collapse tandem duplicates using tandem array IDs from combBed, "
            "choose array representatives using isArrayRep == TRUE where available, "
            "rewrite the filtered TSV, and regenerate the OG list."
        )
    )
    parser.add_argument("--infile", required=True, help="Input PASS pangenes TSV")
    parser.add_argument("--combBed", required=True, help="Input combBed.txt TSV")
    parser.add_argument("--genomes-tsv", required=True, help="Input genomes TSV")
    parser.add_argument("--outfile_filtered", required=True, help="Filtered output TSV")
    parser.add_argument("--outfile_tandems", required=True, help="Tandem report TSV")
    parser.add_argument("--outfile_og_list", required=True, help="Rewritten OG list")
    parser.add_argument(
        "--require-outgroup",
        action="store_true",
        help="Require each retained OG to include at least one outgroup genome"
    )
    args = parser.parse_args()

    infile = Path(args.infile)
    combBed_file = Path(args.combBed)
    genomes_tsv = Path(args.genomes_tsv)
    outfile_filtered = Path(args.outfile_filtered)
    outfile_tandems = Path(args.outfile_tandems)
    outfile_og_list = Path(args.outfile_og_list)

    df = pd.read_csv(infile, sep="\t", dtype=str)
    combBed = pd.read_csv(combBed_file, sep="\t", dtype=str)

    required_pass_cols = ["og", "genome", "chr", "id", "ord"]
    missing_pass = [c for c in required_pass_cols if c not in df.columns]
    if missing_pass:
        raise SystemExit(
            f"Missing required columns in input PASS TSV: {', '.join(missing_pass)}"
        )

    required_combBed_cols = ["id", "genome", "arrayID", "isArrayRep"]
    missing_combBed = [c for c in required_combBed_cols if c not in combBed.columns]
    if missing_combBed:
        raise SystemExit(
            f"Missing required columns in combBed TSV: {', '.join(missing_combBed)}"
        )

    for col in required_pass_cols:
        df[col] = df[col].fillna("").astype(str).str.strip()

    for col in required_combBed_cols:
        combBed[col] = combBed[col].fillna("").astype(str).str.strip()

    df["ord_num"] = pd.to_numeric(df["ord"], errors="coerce")

    input_rows = len(df)
    input_ogs = df["og"].nunique()

    combBed_subset = combBed[["id", "genome", "arrayID", "isArrayRep"]].drop_duplicates()

    duplicated_keys = combBed_subset.duplicated(subset=["id", "genome"], keep=False)
    if duplicated_keys.any():
        dup_df = combBed_subset.loc[duplicated_keys, ["id", "genome"]].drop_duplicates()
        example_rows = dup_df.head(10).to_dict(orient="records")
        raise SystemExit(
            "Duplicate id+genome keys found in combBed, cannot merge unambiguously. "
            f"Example duplicate keys: {example_rows}"
        )

    df = df.merge(
        combBed_subset,
        on=["id", "genome"],
        how="left",
        validate="many_to_one"
    )

    missing_array_info = df["arrayID"].isna().sum()
    missing_rep_info = df["isArrayRep"].isna().sum()

    df["arrayID"] = df["arrayID"].fillna("").astype(str).str.strip()
    df["isArrayRep"] = df["isArrayRep"].fillna("").astype(str).str.strip()
    df["isArrayRep_norm"] = df["isArrayRep"].map(normalize_bool_str)

    kept_rows = []
    tandem_records = []

    n_clusters_using_true_rep = 0
    n_clusters_multiple_true = 0
    n_clusters_fallback = 0

    for (og, genome, chr_id), subdf in df.groupby(["og", "genome", "chr"], sort=False):
        subdf = subdf.sort_values(["ord_num", "id"], na_position="last").copy()

        subdf["collapse_key"] = subdf["arrayID"]
        no_array_mask = subdf["collapse_key"] == ""
        subdf.loc[no_array_mask, "collapse_key"] = (
            "__singleton__" + subdf.loc[no_array_mask, "id"]
        )

        for collapse_key, cluster in subdf.groupby("collapse_key", sort=False):
            cluster = cluster.sort_values(["ord_num", "id"], na_position="last").copy()

            kept, selection_method, n_true_reps = choose_kept_row(cluster)
            kept_rows.append(kept)

            if selection_method == "isArrayRep_TRUE":
                n_clusters_using_true_rep += 1
            elif selection_method == "multiple_TRUE_lowest_ord":
                n_clusters_multiple_true += 1
            elif selection_method == "fallback_lowest_ord":
                n_clusters_fallback += 1

            is_real_array = (
                pd.notna(kept["arrayID"]) and
                str(kept["arrayID"]).strip() != ""
            )

            if is_real_array and len(cluster) > 1:
                removed = cluster.loc[cluster["id"] != kept["id"], "id"].tolist()

                tandem_records.append({
                    "og": og,
                    "genome": genome,
                    "chr": chr_id,
                    "arrayID": kept["arrayID"],
                    "kept_gene": kept["id"],
                    "kept_isArrayRep": kept["isArrayRep"],
                    "selection_method": selection_method,
                    "n_true_reps_in_cluster": n_true_reps,
                    "removed_genes": ",".join(removed),
                    "cluster_size": len(cluster),
                    "n_removed": len(cluster) - 1
                })

    collapsed_df = pd.DataFrame(kept_rows).drop(
        columns=["ord_num", "collapse_key", "isArrayRep_norm"],
        errors="ignore"
    )
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
            columns=[
                "og", "genome", "chr", "arrayID",
                "kept_gene", "kept_isArrayRep", "selection_method",
                "n_true_reps_in_cluster", "removed_genes",
                "cluster_size", "n_removed"
            ]
        ).to_csv(outfile_tandems, sep="\t", index=False)
    else:
        tandem_df.to_csv(outfile_tandems, sep="\t", index=False)

    print("--------------------------------------------------")
    print("collapse_tandems summary")
    print("--------------------------------------------------")
    print(f"Input PASS file: {infile}")
    print(f"Input combBed file: {combBed_file}")
    print(f"Input rows: {input_rows}")
    print(f"Input OGs: {input_ogs}")
    print(f"Rows lacking arrayID after merge: {missing_array_info}")
    print(f"Rows lacking isArrayRep after merge: {missing_rep_info}")
    print(f"Outgroup requirement enabled: {'yes' if args.require_outgroup else 'no'}")
    if args.require_outgroup:
        print(f"Outgroup genomes loaded: {len(outgroup_genomes)}")
    print("--------------------------------------------------")
    print(f"Rows after tandem collapse: {rows_after_collapse}")
    print(f"OGs after tandem collapse: {ogs_after_collapse}")
    print(f"Tandem clusters found: {tandem_clusters_found}")
    print(f"Rows removed by tandem collapse: {rows_removed_by_tandem}")
    print("--------------------------------------------------")
    print(f"Clusters kept via single TRUE array rep: {n_clusters_using_true_rep}")
    print(f"Clusters with multiple TRUE reps, resolved by ord: {n_clusters_multiple_true}")
    print(f"Clusters with no TRUE rep, fell back to ord: {n_clusters_fallback}")
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
