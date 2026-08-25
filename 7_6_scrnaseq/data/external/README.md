# External derived inputs

Small (<2 MB) derived tables produced by analyses outside this repo's scripts,
vendored here so the pipeline is self-contained. Neither is raw sequencing
data — see the main README for that.

- **infercnv_chr8q_scores.csv** — per-cell chr8q CNV scores, from a separate
  inferCNV analysis. Columns: cell barcode (index), `Cell_type`,
  `chr8q_cnv_score`, `chr8q_status`. Consumed by `scripts/04_8q_cnv_score.py`.
- **hM1_8q_blocks.csv** — chr8q gene → block assignments, from the
  `organoid_bulk` analysis. Columns: `symbol`, `block`. Consumed by
  `scripts/04b_8q_block_trajectory.py`.
