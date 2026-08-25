#!/bin/bash
# Downloads ERK-inhibitor time-course RNA-seq (ENA BioProject PRJEB25806) via
# the ENA file-report API. Run from the 7_22_simyc/ repo root.
set -e

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB25806&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
  | tail -n +2 | cut -f2 | tr ';' '\n' | sed '/^$/d' > ena_urls_erk.txt

echo "URL count: $(wc -l < ena_urls_erk.txt)"

mkdir -p fastq_erk
xargs -n 1 -P 8 -I{} wget -q -c https://{} -P fastq_erk/ < ena_urls_erk.txt

echo "Files downloaded: $(ls fastq_erk/ | wc -l)"
