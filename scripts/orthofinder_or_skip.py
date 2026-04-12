#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import sys


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptide-dir", required=True)
    ap.add_argument("--orthofinder-dir", required=True)
    ap.add_argument("--orthofinder-bin", required=True)
    ap.add_argument("--threads", type=int, required=True)
    ap.add_argument("--analysis-threads", type=int, default=1)
    ap.add_argument("--genomes", nargs="+", required=True)
    ap.add_argument("--force", default="false")
    return ap.parse_args()


def as_bool(x: str) -> bool:
    return str(x).lower() in {"1", "true", "yes", "y"}


def read_species_set(of_dir: Path):
    marker = of_dir / "species_set.txt"
    if not marker.exists():
        return None
    return [l.strip() for l in marker.read_text().splitlines() if l.strip()]


def write_species_set(of_dir: Path, genomes):
    marker = of_dir / "species_set.txt"
    marker.write_text("\n".join(sorted(genomes)) + "\n")


def find_results_dir(of_dir: Path):
    """
    OrthoFinder writes:
      <of_dir>/Results_<DATE>/
    We must return the full path to that Results_* directory.
    """
    results = sorted([p for p in of_dir.glob("Results_*") if p.is_dir()])
    if not results:
        return None
    return results[-1]


def write_results_dir_txt(of_dir: Path, results_dir: Path):
    out = of_dir / "results_dir.txt"
    out.write_text(str(results_dir.resolve()) + "\n")


def require_results_dir(of_dir: Path):
    results_dir = find_results_dir(of_dir)
    if results_dir is None:
        raise SystemExit(f"ERROR: could not find Results_* directory inside {of_dir}")
    write_results_dir_txt(of_dir, results_dir)
    return results_dir


def main():
    args = parse_args()

    peptide_dir = Path(args.peptide_dir)
    of_dir = Path(args.orthofinder_dir)
    force = as_bool(args.force)

    if not peptide_dir.exists():
        print(f"ERROR: peptide dir not found: {peptide_dir}", file=sys.stderr)
        sys.exit(2)

    if args.analysis_threads < 1:
        print("ERROR: --analysis-threads must be >= 1", file=sys.stderr)
        sys.exit(2)

    expected = sorted(args.genomes)

    existing = read_species_set(of_dir)
    if of_dir.exists() and existing is not None and sorted(existing) == expected and not force:
        print(f"Found existing OrthoFinder run in {of_dir} with matching species_set.txt; skipping.")
        results_dir = require_results_dir(of_dir)
        print(f"Using existing Results dir: {results_dir}")
        (of_dir / "orthofinder.skipped").write_text("skipped\n")
        return

    print("Running OrthoFinder (new run or mismatched species set).")
    of_dir.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        args.orthofinder_bin,
        "-f", str(peptide_dir),
        "-t", str(args.threads),
        "-a", str(args.analysis_threads),
        "-o", str(of_dir),
    ]

    print("Command:", " ".join(cmd))
    subprocess.check_call(cmd)

    write_species_set(of_dir, expected)

    results_dir = require_results_dir(of_dir)
    print(f"OrthoFinder finished. Results dir: {results_dir}")
    (of_dir / "orthofinder.ran").write_text("ran\n")


if __name__ == "__main__":
    main()
