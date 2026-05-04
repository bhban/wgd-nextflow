#!/bin/bash

# Script to check base directories before downloading and renaming NCBI annotation files
# bash check_ncbi_ftp_dirs.sh ncbi_ftp_base_urls.txt
# ncbi_ftp_base_urls.txt should have one url per line 
# e.g.,
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/939/355/GCA_036939355.2_TProf2.0/
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/187/285/GCA_013187285.2_UArk_CyCryp_3/
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/940/135/GCA_036940135.2_UArk_CyAtom_2/

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 ncbi_ftp_base_urls.txt" >&2
    exit 1
fi

url_file="$1"

if [ ! -s "$url_file" ]; then
    echo "Error: URL file does not exist or is empty: $url_file" >&2
    exit 1
fi

while IFS= read -r url; do
    # Skip blank lines and comments
    if [ -z "$url" ] || [[ "$url" =~ ^[[:space:]]*# ]]; then
        continue
    fi

    # Ensure trailing slash
    url="${url%/}/"

    echo "Checking: $url"

    if wget --spider -q "$url"; then
        echo "  OK: directory exists"
    else
        echo "  FAIL: directory not found or inaccessible"
        echo
        continue
    fi

    page="$(wget -qO - "$url")"

    for suffix in \
        "genomic.gff.gz" \
        "protein.faa.gz" \
        "cds_from_genomic.fna.gz"
    do
        match="$(printf '%s\n' "$page" | grep -oE "GC[AF]_[^\"<> ]+_${suffix}" | head -n 1 || true)"

        if [ -n "$match" ]; then
            echo "  found: $match"
        else
            echo "  missing: *_${suffix}"
        fi
    done

    echo
done < "$url_file"
