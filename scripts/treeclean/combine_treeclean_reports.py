#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine tree-cleaning reports into one summary TSV."
    )
    parser.add_argument("--decomposition-reports", nargs="+", required=True)
    parser.add_argument("--subtree-filter-reports", nargs="+", required=True)
    parser.add_argument("--treeshrink-pruning-report", required=True)
    parser.add_argument("--post-treeshrink-filter-report", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def read_tsv(path: str) -> List[dict]:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_many(paths: Iterable[str]) -> List[dict]:
    rows = []
    for path in paths:
        rows.extend(read_tsv(path))
    return rows


def index_by_unit(rows: List[dict]) -> Dict[str, dict]:
    out = {}
    for row in rows:
        unit = row.get("cleaning_unit_id")
        if unit:
            out[unit] = row
    return out


def decomposition_by_og(rows: List[dict]) -> Dict[str, dict]:
    out = {}
    for row in rows:
        og = row.get("original_og")
        if og and og not in out:
            out[og] = row
    return out


def main() -> None:
    args = parse_args()

    decomp_rows = read_many(args.decomposition_reports)
    subtree_filter_rows = read_many(args.subtree_filter_reports)
    pruning_rows = read_tsv(args.treeshrink_pruning_report)
    post_filter_rows = read_tsv(args.post_treeshrink_filter_report)

    decomp_by_og = decomposition_by_og(decomp_rows)
    subtree_filter_by_unit = index_by_unit(subtree_filter_rows)
    pruning_by_unit = index_by_unit(pruning_rows)
    post_filter_by_unit = index_by_unit(post_filter_rows)

    all_units = sorted(
        set(subtree_filter_by_unit)
        | set(pruning_by_unit)
        | set(post_filter_by_unit)
    )

    fieldnames = [
        "original_og",
        "cleaning_unit_id",
        "subtree_index",
        "n_tips_original_og",
        "n_subtrees_from_original_og",
        "n_cuts_original_og",
        "n_tips_after_decompose",
        "n_species_after_decompose",
        "has_outgroup_after_decompose",
        "passes_filter_after_decompose",
        "fail_reason_after_decompose",
        "n_tips_before_treeshrink",
        "n_removed_by_treeshrink",
        "removed_tips",
        "n_tips_after_treeshrink",
        "n_species_after_treeshrink",
        "has_outgroup_after_treeshrink",
        "passes_filter_after_treeshrink",
        "fail_reason_after_treeshrink",
    ]

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for unit in all_units:
            pre = subtree_filter_by_unit.get(unit, {})
            prune = pruning_by_unit.get(unit, {})
            post = post_filter_by_unit.get(unit, {})

            original_og = (
                pre.get("original_og")
                or prune.get("original_og")
                or post.get("original_og")
                or "NA"
            )

            decomp = decomp_by_og.get(original_og, {})

            writer.writerow(
                {
                    "original_og": original_og,
                    "cleaning_unit_id": unit,
                    "subtree_index": (
                        pre.get("subtree_index")
                        or prune.get("subtree_index")
                        or post.get("subtree_index")
                        or "NA"
                    ),
                    "n_tips_original_og": decomp.get("n_tips_original", "NA"),
                    "n_subtrees_from_original_og": decomp.get("n_subtrees", "NA"),
                    "n_cuts_original_og": decomp.get("n_cuts", "NA"),
                    "n_tips_after_decompose": pre.get("n_tips", "NA"),
                    "n_species_after_decompose": pre.get("n_species", "NA"),
                    "has_outgroup_after_decompose": pre.get("has_outgroup", "NA"),
                    "passes_filter_after_decompose": pre.get("passes_filter", "NA"),
                    "fail_reason_after_decompose": pre.get("fail_reason", "NA"),
                    "n_tips_before_treeshrink": prune.get("n_tips_before", "NA"),
                    "n_removed_by_treeshrink": prune.get("n_removed", "NA"),
                    "removed_tips": prune.get("removed_tips", "NA"),
                    "n_tips_after_treeshrink": post.get("n_tips", prune.get("n_tips_after", "NA")),
                    "n_species_after_treeshrink": post.get("n_species", "NA"),
                    "has_outgroup_after_treeshrink": post.get("has_outgroup", "NA"),
                    "passes_filter_after_treeshrink": post.get("passes_filter", "NA"),
                    "fail_reason_after_treeshrink": post.get("fail_reason", "NA"),
                }
            )


if __name__ == "__main__":
    main()
