#!/usr/bin/env python

"""primary_transcripts.py: Adapted from the orthofinder script. Adds option to pull primary_transcripts for Phytozome genomes."""

import os
import re
import sys
import argparse
from collections import Counter, defaultdict


def CheckFile(fn):
    """
    Checks for:
    - Duplicated accession lines
    """
    accs = set()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                a = l.rstrip()[1:]
                if a in accs:
                    print("\nERROR: duplicated sequence accession:\n%s" % a, file=sys.stderr)
                    print("\nPlease correct this and then rerun the script.\n", file=sys.stderr)
                    return False
                accs.add(a)
    return True


def ScanTags_with_fn(fn, gene_name_fn):
    genes = []
    with open(fn, 'r') as infile:
        for line in infile:
            if not line.startswith(">"):
                continue
            genes.append(gene_name_fn(line))
    print("%d sequences, %d genes" % (len(genes), len(set(genes))), file=sys.stderr)


def GetGeneName_Ensembl(acc_line):
    tokens = [(t.split("=") if "=" in t else t.split(":"))[1]
              for t in acc_line.rstrip().split()
              if ("gene:" in t or "gene=" in t or "locus:" in t or "locus=" in t)]
    if len(tokens) != 1:
        return None
    return tokens[0]


def IsNCBI(fn):
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith(">"):
                l = l.rstrip()
                if l.startswith(">NP_") and l.endswith("]"):
                    return True
                elif l.startswith(">XP_") and l.endswith("]"):
                    return True
                elif l.startswith(">YP_") and l.endswith("]"):
                    return True
                elif l.startswith(">WP_") and l.endswith("]"):
                    return True
                return False
    return False


def GetGeneName_NCBI(acc_line):
    acc_line = acc_line[1:]
    original = acc_line
    acc_line = re.sub("isoform [0-9, A-Z]+ ", "", acc_line)
    acc_line = re.sub("isoform X[0-9, A-Z]+ ", "", acc_line)
    if original != acc_line:
        acc_line = acc_line.split(None, 1)[-1]
    return acc_line


def fasta_first_token(acc_line):
    # acc_line starts with ">"
    return acc_line[1:].rstrip().split(None, 1)[0]


def parse_gff_attr(attr_field):
    d = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def phytozome_keep_ids_from_gff(gff_fn):
    """
    Return a set of IDs for which longest=1 (keep),
    and longest=0 (drop) is ignored.
    """
    keep = set()
    n_seen = 0
    n_keep = 0
    with open(gff_fn, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in ("mRNA", "transcript"):
                continue

            attrs = parse_gff_attr(parts[8])
            if "longest" not in attrs:
                continue

            n_seen += 1
            longest_val = attrs.get("longest", "").strip()
            if longest_val != "1":
                continue

            # Prefer ID=..., else Parent=..., else Name=...
            tid = attrs.get("ID") or attrs.get("Name") or attrs.get("Parent")
            if not tid:
                continue

            keep.add(tid)
            n_keep += 1

    print(f"Phytozome GFF: saw {n_seen} transcript features with longest=; keeping {n_keep} with longest=0", file=sys.stderr)
    return keep


def CreatePrimaryTranscriptsFile(fn, dout, gene_name_fn, q_use_original_accession_line):
    max_gene_lens = defaultdict(int)
    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1
    nAcc = 0
    nGeneUnidentified = 0
    acc_to_use = defaultdict(str)
    iLine = -1
    while iLine < N:
        iLine += 1
        line = lines[iLine]
        if not line.startswith(">"):
            continue
        nAcc += 1
        iLineAcc = iLine
        gene = gene_name_fn(line)
        if gene is None:
            nGeneUnidentified += 1
            continue
        l = 0
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if line.startswith(">"):
                iLine -= 1
                break
            l += len(line.rstrip())
        if l > max_gene_lens[gene]:
            max_gene_lens[gene] = l
            acc_to_use[gene] = iLineAcc

    print("Found %d accessions, %d genes, %d unidentified transcripts" % (nAcc, len(max_gene_lens), nGeneUnidentified), file=sys.stderr)

    nGenesWriten = 0
    outfn = os.path.join(dout, os.path.basename(fn))
    with open(outfn, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"):
                continue
            gene = gene_name_fn(line)
            if gene is not None and iLine != acc_to_use[gene]:
                continue
            if q_use_original_accession_line or gene is None:
                acc_line_out = line + "\n"
            else:
                acc_line_out = ">%s\n" % gene
            nGenesWriten += 1
            outfile.write(acc_line_out)
            while iLine < N:
                iLine += 1
                line = lines[iLine]
                if line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(line + "\n")

    print("Wrote %d genes" % nGenesWriten, file=sys.stderr)
    if nGenesWriten != len(max_gene_lens) + nGeneUnidentified:
        print("ERROR", file=sys.stderr)
        raise Exception
    print(outfn, file=sys.stderr)


def CreatePrimaryTranscriptsFile_Phytozome(fn, dout, keep_ids):
    """
    Keep isoforms whose FASTA first token matches a transcript ID with longest=0 in the GFF.
    Preserves original full header lines.
    """
    with open(fn, 'r') as infile:
        lines = [l.rstrip() for l in infile]
    N = len(lines) - 1

    nAcc = 0
    nKept = 0
    nDropped = 0
    nNoMatch = 0

    outfn = os.path.join(dout, os.path.basename(fn))
    with open(outfn, 'w') as outfile:
        iLine = -1
        while iLine < N:
            iLine += 1
            line = lines[iLine]
            if not line.startswith(">"):
                continue

            nAcc += 1
            tok = fasta_first_token(line)

            if tok not in keep_ids:
                # Track whether it's unknown vs explicitly not kept
                nDropped += 1 if tok in keep_ids or tok else 1
                if tok not in keep_ids:
                    nNoMatch += 1
                # skip sequence body
                while iLine < N:
                    iLine += 1
                    if lines[iLine].startswith(">"):
                        iLine -= 1
                        break
                continue

            nKept += 1
            outfile.write(line + "\n")
            while iLine < N:
                iLine += 1
                seq_line = lines[iLine]
                if seq_line.startswith(">"):
                    iLine -= 1
                    break
                outfile.write(seq_line + "\n")

    print(f"Found {nAcc} accessions; kept {nKept}; skipped {nAcc - nKept}", file=sys.stderr)
    print(f"Output: {outfn}", file=sys.stderr)


def last_dot(text):
    return text[1:].rstrip().rsplit(".", 1)[0]


def space(text):
    return text[1:].rstrip().split(None, 1)[0]


function_dict = {"last_dot": last_dot, "space": space}


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    ap = argparse.ArgumentParser(
        description="Create primary transcript FASTA. Default auto-detects NCBI vs Ensembl behavior. "
                    "Phytozome mode uses GFF longest= flag."
    )
    ap.add_argument("pep_fasta", help="Protein FASTA input")
    ap.add_argument("--mode", choices=["auto", "ensembl", "ncbi", "phytozome"], default="auto")
    ap.add_argument("--phytozome-gff", default=None, help="Required if --mode phytozome; GFF3 file to read longest=")
    ap.add_argument("--gene-fn", choices=["last_dot", "space"], default=None,
                    help="Optional override for gene grouping function (legacy behavior)")
    args = ap.parse_args(argv)

    fn = args.pep_fasta

    if not CheckFile(fn):
        return

    dout = os.path.join(os.path.dirname(os.path.abspath(fn)), "primary_transcripts")
    if not os.path.exists(dout):
        os.mkdir(dout)

    # Legacy override: allow explicit grouping fn
    if args.gene_fn is not None:
        gene_name_function = function_dict[args.gene_fn]
        ScanTags_with_fn(fn, gene_name_function)
        CreatePrimaryTranscriptsFile(fn, dout, gene_name_function, q_use_original_accession_line=False)
        return

    if args.mode == "phytozome":
        if not args.phytozome_gff:
            raise SystemExit("ERROR: --phytozome-gff is required when --mode phytozome")
        keep_ids = phytozome_keep_ids_from_gff(args.phytozome_gff)
        CreatePrimaryTranscriptsFile_Phytozome(fn, dout, keep_ids)
        return

    if args.mode == "ncbi" or (args.mode == "auto" and IsNCBI(fn)):
        print("Identified as NCBI file", file=sys.stderr)
        gene_name_function = GetGeneName_NCBI
        q_use_original_accession_line = True
    else:
        gene_name_function = GetGeneName_Ensembl
        q_use_original_accession_line = False
        print('Looking for "gene=" or "gene:" to identify isoforms of same gene', file=sys.stderr)

    CreatePrimaryTranscriptsFile(fn, dout, gene_name_function, q_use_original_accession_line)


if __name__ == "__main__":
    main()
