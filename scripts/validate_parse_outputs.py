#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genespace-wd", required=True, help="e.g. genespace/workingDir")
    ap.add_argument("--genomes", nargs="+", required=True, help="genome IDs")
    return ap.parse_args()


def main():
    args = parse_args()
    wd = Path(args.genespace_wd)

    bed_dir = wd / "bed"
    pep_dir = wd / "peptide"

    missing = []
    for g in args.genomes:
        bed = bed_dir / f"{g}.bed"
        fa = pep_dir / f"{g}.fa"
        if not bed.exists():
            missing.append(str(bed))
        if not fa.exists():
            missing.append(str(fa))

    if missing:
        print("ERROR: parse_annotations outputs missing:")
        for p in missing:
            print(p)
        sys.exit(2)

    print(f"OK: found {len(args.genomes)} bed files in {bed_dir}")
    print(f"OK: found {len(args.genomes)} fa files in {pep_dir}")


if __name__ == "__main__":
    main()
