#!/bin/bash
set -euo pipefail

manifest="${1:-annotation_urls.tsv}"

mkdir -p gff protein cds

tail -n +2 "$manifest" | while IFS=$'\t' read -r genome_id file_type url outfile; do
    echo "Downloading ${outfile}"

    case "$file_type" in
        gff3)
            wget -qO - "$url" | gunzip -c > "gff/${outfile}"
            ;;
        pep)
            wget -qO - "$url" | gunzip -c > "protein/${outfile}"
            ;;
        cds)
            wget -qO - "$url" | gunzip -c > "cds/${outfile}"
            ;;
        *)
            echo "Skipping unknown file_type: ${file_type}" >&2
            ;;
    esac
done
