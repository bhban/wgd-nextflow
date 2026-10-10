#!/usr/bin/env python3
"""Build uncollapsed protein/CDS OG FASTAs from OrthoFinder Orthogroups.tsv.

Strictly matches gene IDs across OrthoFinder and raw per-species FASTAs;
no GENESPACE membership or synteny information is consulted.
"""
import argparse
import csv
import gzip
import re
from pathlib import Path


def records(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as fh:
        name, chunks = None, []
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
        if name is not None:
            yield name, ''.join(chunks)


def index_fasta(path):
    d = {}
    for key, seq in records(path):
        if key in d:
            raise ValueError(f'Duplicate FASTA identifier {key} in {path}')
        d[key] = seq
    return d


def locate(directory, species, suffix):
    matches = sorted(Path(directory).glob(f'{species}.{suffix}'))
    if not matches:
        raise FileNotFoundError(f'No {species}.{suffix} in {directory}')
    return matches[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--orthofinder-dir', required=True)
    p.add_argument('--genomes-tsv', required=True)
    p.add_argument('--protein-dir', required=True)
    p.add_argument('--cds-dir', required=True)
    p.add_argument('--protein-ext', required=True)
    p.add_argument('--cds-ext', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--min-species', type=int, default=4)
    a = p.parse_args()

    root = Path(a.orthofinder_dir)
    found = sorted(root.rglob('Orthogroups.tsv'))
    if not found:
        raise FileNotFoundError(f'Orthogroups.tsv not found under {root}')
    if len(found) > 1:
        raise RuntimeError(f'Multiple Orthogroups.tsv files found: {found}; pass a single OrthoFinder results directory')
    og_table = found[0]
    with open(a.genomes_tsv, newline='') as fh:
        genomes = [x['genome_id'].strip() for x in csv.DictReader(fh, delimiter='\t')]
    out = Path(a.outdir)
    (out / 'og_fasta_aa').mkdir(parents=True, exist_ok=True)
    (out / 'og_fasta_cds').mkdir(parents=True, exist_ok=True)
    seqs = {}
    for sp in genomes:
        seqs[sp] = (
            index_fasta(locate(a.protein_dir, sp, a.protein_ext)),
            index_fasta(locate(a.cds_dir, sp, a.cds_ext)),
        )
    n_written = 0
    with open(og_table, newline='') as fh, (out / 'og_report.tsv').open('w') as report:
        rd = csv.DictReader(fh, delimiter='\t')
        if not rd.fieldnames or rd.fieldnames[0] != 'Orthogroup':
            raise ValueError('Expected first Orthogroups.tsv column Orthogroup')
        cols = rd.fieldnames[1:]
        if set(cols) != set(genomes):
            raise ValueError(f'OrthoFinder species columns do not match genomes TSV. OrthoFinder={cols}; genomes={genomes}')
        report.write('orthogroup\tspecies_count\tgene_count\tstatus\n')
        for row in rd:
            og = row['Orthogroup']
            genes = [(sp, gene.strip()) for sp in cols for gene in (row[sp] or '').split(',') if gene.strip()]
            represented = len(set(sp for sp, _ in genes))
            status = 'OK' if represented >= a.min_species else 'TOO_FEW_SPECIES'
            aa, nt = [], []
            if status == 'OK':
                for sp, gene in genes:
                    aa_dict, nt_dict = seqs[sp]
                    if gene not in aa_dict or gene not in nt_dict:
                        status = 'MISSING_SEQUENCE'
                        break
                    # Gene identifiers are preserved; species prefix avoids ambiguous tips.
                    label = f'{sp}|{gene}'
                    aa.append(f'>{label}\n{aa_dict[gene]}\n')
                    nt.append(f'>{label}\n{nt_dict[gene]}\n')
            if status == 'OK':
                safe = re.sub(r'[^A-Za-z0-9_.-]', '_', og)
                (out / 'og_fasta_aa' / f'{safe}.fasta').write_text(''.join(aa))
                (out / 'og_fasta_cds' / f'{safe}.fasta').write_text(''.join(nt))
                n_written += 1
            report.write(f'{og}\t{represented}\t{len(genes)}\t{status}\n')
    if n_written == 0:
        raise RuntimeError('No orthogroups written; inspect og_report.tsv and identifier mappings')
    print(f'Wrote {n_written} raw OrthoFinder orthogroups from {og_table}')


if __name__ == '__main__':
    main()
