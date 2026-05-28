#!/usr/bin/env python3

import argparse
import os
import glob
import pandas as pd
from ete3 import Tree


def parse_support(node):
    """
    Return branch support as a numeric value.

    Priority:
    1. Use numeric internal node labels, e.g. )87:0.02
    2. Use node.support only if it appears to be meaningful
    3. Return None for missing support

    Important:
    ETE can assign node.support = 1.0 by default when no support
    value is present or when support labels are parsed as node names.
    Therefore, unlabelled nodes with support == 1.0 are treated as missing.
    """

    if node.name not in [None, ""]:
        try:
            return float(node.name)
        except ValueError:
            return None

    if node.support is not None:
        try:
            support = float(node.support)
        except ValueError:
            return None

        if support == 1.0:
            return None

        return support

    return None


def get_support_values(tree_file, tree_format):
    """
    Extract support values from all internal nodes in a tree.
    """

    tree = Tree(tree_file, format=tree_format)

    rows = []

    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue

        support = parse_support(node)
        descendant_tips = node.get_leaf_names()

        rows.append({
            "tree_file": os.path.basename(tree_file),
            "node_name": node.name if node.name else "",
            "support": support,
            "n_descendant_tips": len(descendant_tips),
            "descendant_tips": ",".join(descendant_tips)
        })

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Summarise branch support values across rooted Newick trees."
    )

    parser.add_argument(
        "--rooted-dir",
        required=True,
        help="Directory containing rooted Newick trees"
    )

    parser.add_argument(
        "--pattern",
        default="*.nwk",
        help="File pattern for rooted trees. Default: *.nwk"
    )

    parser.add_argument(
        "--tree-format",
        type=int,
        default=1,
        help="ETE tree format. Default: 1"
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix"
    )

    args = parser.parse_args()

    tree_files = sorted(glob.glob(os.path.join(args.rooted_dir, args.pattern)))

    if len(tree_files) == 0:
        raise FileNotFoundError(
            f"No trees found in {args.rooted_dir} matching pattern {args.pattern}"
        )

    all_rows = []

    for tree_file in tree_files:
        try:
            rows = get_support_values(tree_file, args.tree_format)
            all_rows.extend(rows)
        except Exception as e:
            all_rows.append({
                "tree_file": os.path.basename(tree_file),
                "node_name": "",
                "support": None,
                "n_descendant_tips": None,
                "descendant_tips": "",
                "error": str(e)
            })

    df = pd.DataFrame(all_rows)

    if "error" not in df.columns:
        df["error"] = ""

    node_out = f"{args.out_prefix}.node_support.tsv"
    df.to_csv(node_out, sep="\t", index=False)

    clean = df[df["error"].fillna("") == ""].copy()
    clean = clean.dropna(subset=["support"])

    if clean.empty:
        print(f"Read {len(tree_files)} trees")
        print(f"Wrote {node_out}")
        print("No labelled branch support values were found.")
        return

    clean["support_bin"] = pd.cut(
        clean["support"],
        bins=[-0.01, 50, 70, 90, 95, 100],
        labels=["<=50", "50-70", "70-90", "90-95", "95-100"],
        include_lowest=True
    )

    per_tree = (
        clean
        .groupby("tree_file")
        .agg(
            n_labelled_internal_nodes=("support", "count"),
            mean_support=("support", "mean"),
            median_support=("support", "median"),
            min_support=("support", "min"),
            max_support=("support", "max"),
            n_support_lt_50=("support", lambda x: (x < 50).sum()),
            n_support_lt_70=("support", lambda x: (x < 70).sum()),
            n_support_lt_90=("support", lambda x: (x < 90).sum()),
            n_support_ge_95=("support", lambda x: (x >= 95).sum())
        )
        .reset_index()
    )

    per_tree["prop_support_lt_50"] = (
        per_tree["n_support_lt_50"] / per_tree["n_labelled_internal_nodes"]
    )
    per_tree["prop_support_lt_70"] = (
        per_tree["n_support_lt_70"] / per_tree["n_labelled_internal_nodes"]
    )
    per_tree["prop_support_lt_90"] = (
        per_tree["n_support_lt_90"] / per_tree["n_labelled_internal_nodes"]
    )
    per_tree["prop_support_ge_95"] = (
        per_tree["n_support_ge_95"] / per_tree["n_labelled_internal_nodes"]
    )

    tree_out = f"{args.out_prefix}.per_tree_support_summary.tsv"
    per_tree.to_csv(tree_out, sep="\t", index=False)

    missing_by_tree = (
        df[df["error"].fillna("") == ""]
        .groupby("tree_file")
        .agg(
            n_internal_nodes_total=("tree_file", "count"),
            n_missing_support=("support", lambda x: x.isna().sum())
        )
        .reset_index()
    )

    missing_by_tree["prop_missing_support"] = (
        missing_by_tree["n_missing_support"] /
        missing_by_tree["n_internal_nodes_total"]
    )

    missing_out = f"{args.out_prefix}.missing_support_summary.tsv"
    missing_by_tree.to_csv(missing_out, sep="\t", index=False)

    overall = clean["support"].describe().to_frame().T
    overall["n_trees"] = clean["tree_file"].nunique()
    overall["n_labelled_internal_nodes"] = clean.shape[0]
    overall["n_internal_nodes_total"] = df[df["error"].fillna("") == ""].shape[0]
    overall["n_missing_support"] = df[
        (df["error"].fillna("") == "") & df["support"].isna()
    ].shape[0]

    overall["prop_missing_support"] = (
        overall["n_missing_support"] / overall["n_internal_nodes_total"]
    )

    overall["n_support_lt_50"] = (clean["support"] < 50).sum()
    overall["n_support_lt_70"] = (clean["support"] < 70).sum()
    overall["n_support_lt_90"] = (clean["support"] < 90).sum()
    overall["n_support_ge_95"] = (clean["support"] >= 95).sum()

    overall["prop_support_lt_50"] = (
        overall["n_support_lt_50"] / overall["n_labelled_internal_nodes"]
    )
    overall["prop_support_lt_70"] = (
        overall["n_support_lt_70"] / overall["n_labelled_internal_nodes"]
    )
    overall["prop_support_lt_90"] = (
        overall["n_support_lt_90"] / overall["n_labelled_internal_nodes"]
    )
    overall["prop_support_ge_95"] = (
        overall["n_support_ge_95"] / overall["n_labelled_internal_nodes"]
    )

    overall_out = f"{args.out_prefix}.overall_support_summary.tsv"
    overall.to_csv(overall_out, sep="\t", index=False)

    bins = (
        clean
        .groupby("support_bin", observed=False)
        .size()
        .reset_index(name="n_nodes")
    )

    bins["prop_nodes"] = bins["n_nodes"] / bins["n_nodes"].sum()

    bins_out = f"{args.out_prefix}.support_bins.tsv"
    bins.to_csv(bins_out, sep="\t", index=False)

    error_rows = df[df["error"].fillna("") != ""].copy()

    if not error_rows.empty:
        error_out = f"{args.out_prefix}.errors.tsv"
        error_rows.to_csv(error_out, sep="\t", index=False)
        print(f"Wrote {error_out}")

    print(f"Read {len(tree_files)} trees")
    print(f"Wrote {node_out}")
    print(f"Wrote {tree_out}")
    print(f"Wrote {missing_out}")
    print(f"Wrote {overall_out}")
    print(f"Wrote {bins_out}")


if __name__ == "__main__":
    main()
