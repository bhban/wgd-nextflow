#!/usr/bin/env python3
import argparse
from collections import Counter


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-gff", required=True)
    ap.add_argument("--out-gff", required=True)
    ap.add_argument("--chr-dict", required=True)
    ap.add_argument("--strict", action="store_true",
                    help="Fail if any seqid in the GFF is not present in the mapping.")
    ap.add_argument("--report-top", type=int, default=20,
                    help="How many unmapped seqids to report (default 20).")
    return ap.parse_args()


def looks_like_int(s: str) -> bool:
    # allow plain digits only for detecting BED coord columns
    return s.isdigit()


def load_chr_dict(tsv_path: str):
    """
    Supports two input styles:

    Style A (2+ columns):
      old_seqid <TAB> new_seqid [<TAB> ...]
      (header optional)

    Style B (chr_lengths.bed / chrom sizes BED-ish):
      #chrom <TAB> chromStart <TAB> chromEnd <TAB> name
      chr1   0   88077321   1
      NC_... 0   9953      NC_...

    For Style B, interpret mapping as:
      old_seqid = column 4 (name)
      new_seqid = column 1 (#chrom)

    No integer-to-chr logic is applied.
    """
    mapping = {}

    with open(tsv_path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            # Skip comment lines (but allow "#chrom" header line to be skipped by header logic below)
            if parts[0].startswith("#") and parts[0] != "#chrom":
                continue

            # Style B header
            if parts[0] == "#chrom":
                continue

            # Style B data row detection: 4+ cols and cols 2/3 look like numeric coordinates
            if len(parts) >= 4 and looks_like_int(parts[1]) and looks_like_int(parts[2]):
                new_seqid = parts[0]
                old_seqid = parts[3]

                if old_seqid and new_seqid:
                    mapping[old_seqid] = new_seqid
                continue

            # Otherwise treat as Style A (old <TAB> new)
            # Skip common headers
            if parts[0].lower() in ("old_seqid", "old", "seqid") and parts[1].lower() in ("new_seqid", "new"):
                continue

            old_seqid, new_seqid = parts[0], parts[1]
            if old_seqid and new_seqid:
                mapping[old_seqid] = new_seqid

    if not mapping:
        raise SystemExit(f"ERROR: no mappings loaded from {tsv_path}")

    return mapping


def main():
    args = parse_args()
    mapping = load_chr_dict(args.chr_dict)

    n_feature_lines = 0
    n_mapped = 0
    n_changed = 0
    unmapped = Counter()

    with open(args.in_gff, "r") as fin, open(args.out_gff, "w") as fout:
        for line in fin:
            if not line or line.startswith("#"):
                fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue

            n_feature_lines += 1
            sid = parts[0]

            if sid in mapping:
                n_mapped += 1
                new_sid = mapping[sid]
                if new_sid != sid:
                    parts[0] = new_sid
                    n_changed += 1
                # if equal, we still count it as mapped but not changed
                fout.write("\t".join(parts) + "\n")
            else:
                unmapped[sid] += 1
                fout.write(line)

    print(f"Loaded {len(mapping)} mappings from: {args.chr_dict}")
    print(f"GFF feature lines: {n_feature_lines}")
    print(f"Seqids matched mapping keys: {n_mapped}/{n_feature_lines}")
    print(f"Seqids actually changed: {n_changed}/{n_feature_lines}")
    print(f"Output: {args.out_gff}")

    if unmapped:
        print(f"Unmapped seqids (top {args.report_top} by count):")
        for sid, cnt in unmapped.most_common(args.report_top):
            print(f"  {sid}\t{cnt}")

        if args.strict:
            raise SystemExit("ERROR: unmapped seqids found and --strict was set.")


if __name__ == "__main__":
    main()
