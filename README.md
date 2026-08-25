# PDAC scRNA-seq: MYC / glycosylation analysis

Analysis of the PDAC (pancreatic ductal adenocarcinoma) single-cell RNA-seq
dataset `StdWf1_PRJCA001063_CRC_besca2`, examining MYC and glycosyltransferase
expression, chr8q amplification signal, and marker-gene panels across cell
types.

## Setup

```bash
pip install -r requirements.txt
```

Tested with Python 3.13.11. `06_heatmaps_8q.py` queries Ensembl BioMart on its
first run (needs network access) and caches the result locally afterward.

## Data

**Raw input** (not included in this repo — too large for git):
`StdWf1_PRJCA001063_CRC_besca2.raw.h5ad`, expected at the repo root. Source
accession: **PRJCA001063** 

**External derived inputs** (vendored in `data/external/`, see
[`data/external/README.md`](data/external/README.md) for provenance):
- `data/external/infercnv_chr8q_scores.csv`
- `data/external/hM1_8q_blocks.csv`

**Generated artifacts** (gitignored, rebuilt by running the pipeline):
`adata_relabeled.h5ad`, `processed_annotated.h5ad`, `cache/*.h5ad`,
`cache/*.pkl`, `results/figures/*.pdf`. `cache/genes_8q_biomart.csv` (the
BioMart chr8q gene list) is small and *is* committed, so `06_heatmaps_8q.py`
doesn't need network access on a fresh checkout.

## Pipeline

Run from the repo root (`7_6_scrnaseq/`), e.g. `python3 scripts/01b_quick_relabel.py`.
Numbering reflects dependency order, not necessarily one-time-only —
`01`/`01b` are the two preprocessing entry points everything else builds on.

| Script | Input | Output | Purpose |
|---|---|---|---|
| `01_qc_clustering_annotation.py` | raw h5ad | `processed_annotated.h5ad` | Full QC, PCA/UMAP, clustering, derived marker genes. Keeps a raw `counts` layer — required by `05`/`05b`. |
| `01b_quick_relabel.py` | raw h5ad | `adata_relabeled.h5ad` | Lightweight normalize + log1p + cell-type relabel, for scripts that only need expression values (no UMAP/clustering). |
| `02_compute_specificity_stats.py` | `adata_relabeled.h5ad` | `cache/adata_t_processed.h5ad`, `cache/pts.pkl` | Tumor-only Wilcoxon rank-genes (slow step), cached for plotting scripts. |
| `03_myc_glyco_panel.py` | `adata_relabeled.h5ad` | figures | MYC/glycotransferase panel: dotplot, violins, two heatmap styles. |
| `03b_glyco_pathway_panel.py` | `adata_relabeled.h5ad` | figures | Core glycosyltransferase panel: violins + dotplot. |
| `03c_nglycan_oglycan_panel.py` | `adata_relabeled.h5ad` | figures | Broader N-/O-glycan pathway gene panel: violins + dotplot. |
| `03d_b3gnt3_focused.py` | `adata_relabeled.h5ad` | figures | B3GNT3 single-gene dotplot, violin, two heatmap styles. |
| `03e_cbioportal_top15.py` | `adata_relabeled.h5ad` | figures | Top-15 cBioPortal-altered genes (+B3GNT3) dotplot. |
| `04_8q_cnv_score.py` | `data/external/infercnv_chr8q_scores.csv` | figures | Malignant vs. immune chr8q CNV score comparison (Mann-Whitney). |
| `04b_8q_block_trajectory.py` | `processed_annotated.h5ad`, `data/external/hM1_8q_blocks.csv` | figures | Per-8q-block dotplots grouped by tumor/normal condition. |
| `05_myc_gt_coexpression.py` | `processed_annotated.h5ad` | figures | MYC/glycotransferase co-detection UpSet plot (malignant cells). |
| `05b_coexpression_sensitivity.py` | `processed_annotated.h5ad` | figures | MYC + top-2-by-Jaccard glycotransferase co-expression vs. detection threshold. |
| `06_heatmaps_8q.py` | `cache/adata_t_processed.h5ad`, `cache/pts.pkl` | figures | Per-cell heatmap of all chr8q genes (BioMart-derived gene list, cached). |


All figures are written to `results/figures/`.

## Notes on reproducibility

- `01_qc_clustering_annotation.py` seeds PCA/neighbors/UMAP (`random_state=0`).
- The gene-symbol alias `TSTA3` → `GFUS` and the cell-type relabeling
  (`Ductal cell type 2` → `Malignant ductal cells`, etc.) are centralized in
  `scripts/_common.py` so every script applies them identically.
- `05b_coexpression_sensitivity.py` pins the two glycotransferases compared
  against MYC (`GFUS`, `GPAA1`, chosen by Jaccard overlap) rather than
  re-deriving them on every run, so the published figure's gene identities
  don't silently change if the input data changes. Set `PIN_GT1_GT2 = False`
  in that script to re-derive them instead.
