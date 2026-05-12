#!/bin/bash
set -euo pipefail

mkdir -p fasta

tail -n +2 assembly_urls.tsv | while IFS=$'\t' read -r genome_id url outfile; do
    echo "Downloading ${outfile}"
    wget -qO "fasta/${outfile}" "${url}"
done
