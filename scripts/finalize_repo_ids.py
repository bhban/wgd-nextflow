#!/usr/bin/env python3
import argparse
import re


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, choices=["phytozome", "ensembl", "ncbi"])
    ap.add_argument("--in-gff", required=True)
    ap.add_argument("--out-gff", required=True)
    ap.add_argument("--in-pep", required=True)
    ap.add_argument("--out-pep", required=True)
    return ap.parse_args()


def strip_final_dot_number(s: str) -> str:
    return re.sub(r"\.\d+$", "", s)


def parse_attr_field(attr: str):
    items = []
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            items.append((k, v))
        else:
            items.append((part, None))
    return items


def rebuild_attr_field(items):
    out = []
    for k, v in items:
        if v is None:
            out.append(k)
        else:
            out.append(f"{k}={v}")
    return ";".join(out)


def rewrite_gff_phytozome(in_gff: str, out_gff: str):
    n_changed = 0
    with open(in_gff, "r") as fin, open(out_gff, "w") as fout:
        for line in fin:
            if not line or line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue

            if parts[2] != "mRNA":
                fout.write(line)
                continue

            attrs = parse_attr_field(parts[8])
            new_attrs = []
            for k, v in attrs:
                if k == "ID" and v is not None:
                    v2 = strip_final_dot_number(v)
                    if v2 != v:
                        n_changed += 1
                    new_attrs.append((k, v2))
                else:
                    new_attrs.append((k, v))

            parts[8] = rebuild_attr_field(new_attrs)
            fout.write("\t".join(parts) + "\n")

    return n_changed


def rewrite_gff_passthrough(in_gff: str, out_gff: str):
    with open(in_gff, "r") as fin, open(out_gff, "w") as fout:
        for line in fin:
            fout.write(line)


def rewrite_pep_strip_final_dot_number(in_pep: str, out_pep: str):
    n_changed = 0
    with open(in_pep, "r") as fin, open(out_pep, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue

            header = line[1:].rstrip("\n")
            if " " in header:
                tok, rest = header.split(" ", 1)
                rest = " " + rest
            else:
                tok, rest = header, ""

            tok2 = strip_final_dot_number(tok)
            if tok2 != tok:
                n_changed += 1

            fout.write(f">{tok2}{rest}\n")

    return n_changed


def rewrite_pep_passthrough(in_pep: str, out_pep: str):
    with open(in_pep, "r") as fin, open(out_pep, "w") as fout:
        for line in fin:
            fout.write(line)


def main():
    args = parse_args()

    # GFF handling
    if args.source == "phytozome":
        n_gff = rewrite_gff_phytozome(args.in_gff, args.out_gff)
    else:
        rewrite_gff_passthrough(args.in_gff, args.out_gff)
        n_gff = 0

    # PEP handling
    if args.source in ("phytozome", "ensembl"):
        n_pep = rewrite_pep_strip_final_dot_number(args.in_pep, args.out_pep)
    else:
        rewrite_pep_passthrough(args.in_pep, args.out_pep)
        n_pep = 0

    print(f"source={args.source}")
    print(f"gff_id_changes={n_gff}")
    print(f"pep_header_changes={n_pep}")


if __name__ == "__main__":
    main()
