#!/bin/bash
# Downloads siMYC knockdown RNA-seq (SRA BioProject PRJNA1018107) via sra-tools.
# Run from the 7_22_simyc/ repo root. Requires: sra-tools (prefetch, fasterq-dump).
set -e

mkdir -p sra_cache fastq

cat data/SRR_list.txt | xargs -n 1 -P 4 -I{} prefetch {} -O sra_cache/
cat data/SRR_list.txt | xargs -n 1 -P 4 -I{} fasterq-dump --split-files -e 4 --outdir fastq/ sra_cache/{}/{}.sra

gzip fastq/*.fastq

echo "Files downloaded: $(ls fastq/ | wc -l)"
