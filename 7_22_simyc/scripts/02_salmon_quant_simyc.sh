#!/bin/bash
# Salmon quantification for the siMYC knockdown samples.
#
# NOTE: no original script for this step was found alongside the repo's other
# download/quant scripts (only the ERK-side equivalent, 02b, existed) — this
# was reconstructed to match that script's pattern. Verify it matches how
# salmon_quant_simyc/ was actually produced before relying on it to
# regenerate published results.
#
# Reference index: built from GENCODE v46 (GRCh38 primary assembly) — not
# included in this repo. Point SALMON_INDEX at your own copy, or set the
# env var.
SALMON_INDEX="${SALMON_INDEX:-../reference/human/salmon_index_human}"
FASTQ_DIR=fastq
OUT_DIR=salmon_quant_simyc
mkdir -p "$OUT_DIR"

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
  SAMPLE=$(basename "$R1" _1.fastq.gz)
  R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"

  if [ -d "$OUT_DIR/$SAMPLE" ]; then
    echo "Skipping $SAMPLE, already done"
    continue
  fi

  salmon quant -i "$SALMON_INDEX" -l A -1 "$R1" -2 "$R2" -p 8 --validateMappings -o "$OUT_DIR/$SAMPLE" \
    || echo "FAILED: $SAMPLE" >> "$OUT_DIR/failed_samples.txt"
done

echo "Quantification complete."
