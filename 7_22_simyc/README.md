# siMYC / ERK-inhibitor bulk RNA-seq

Two bulk RNA-seq datasets compared on a shared panel of MYC/glycosylation
genes (`EXT1`, `ST3GAL1`, `GPAA1`, `HAS2`, `MYC`, `GFUS`):

1. **siMYC** — MYC knockdown vs. control across 8 pancreatic cancer cell lines.
2. **ERKi** — an ERK-inhibitor time course (vehicle, 1h, 4h, 12h, 24h) across 6 cell lines.

## Setup

R packages (tested versions): `tximport` 1.26.1, `DESeq2` 1.38.3,
`ggplot2` 3.5.2, `dplyr` 1.2.1. Also requires `salmon` and `sra-tools`
(`prefetch`, `fasterq-dump`) on PATH for the download/quant steps.

## Data

**Raw input** (not included — too large for git):
- siMYC: SRA BioProject **PRJNA1018107**, run accessions in `data/SRR_list.txt` / `data/runinfo.csv`.
- ERKi: ENA BioProject **PRJEB25806**, run accessions in `data/erk_metadata.csv` / `data/runinfo_erk.csv`.

**External reference** (not included — build/obtain separately):
- Salmon index + tx2gene mapping built from GENCODE v46 (GRCh38 primary
  assembly). Scripts default to `../reference/human/{salmon_index_human,tx2gene.tsv}`
  relative to this folder; override with the `SALMON_INDEX` / `TX2GENE_PATH`
  env vars.
- `cptac_glyco/` (`PDAC_meta.xlsx`, `PDAC_nglycoform.tsv`) — CPTAC PDAC
  proteogenomic dataset, glycoproteomic arm. Cao L, Huang C, Zhou DC, et al.
  "Proteogenomic characterization of pancreatic ductal adenocarcinoma."
  Cell. 2021;184(19):5031-5052.e26. doi:10.1016/j.cell.2021.08.023.
  Processed tables available via LinkedOmics
  (linkedomics.org/data_download/CPTAC-PDAC/) or the Proteomic Data Commons
  (pdc.cancer.gov).

**Generated artifacts** (gitignored, rebuilt by running the pipeline):
`fastq/`, `fastq_erk/`, `sra_cache/`, `salmon_quant_simyc/`,
`salmon_quant_erk/`, `results/`.

## Pipeline

Run from this folder (`7_22_simyc/`).

| Script | Input | Output | Purpose |
|---|---|---|---|
| `01_pull_simyc_fastq.sh` | `data/SRR_list.txt` | `fastq/` | Download siMYC reads via sra-tools. |
| `01b_pull_erk_fastq.sh` | — (queries ENA) | `fastq_erk/` | Download ERKi reads via the ENA API. |
| `02_salmon_quant_simyc.sh` | `fastq/` | `salmon_quant_simyc/` | Salmon quant, siMYC samples. |
| `02b_salmon_quant_erk.sh` | `fastq_erk/` | `salmon_quant_erk/` | Salmon quant, ERKi samples. |
| `03_simyc_deseq2.R` | `salmon_quant_simyc/`, `data/simyc_metadata.csv` | `results/siMYC_DE_results.csv` | DESeq2: siMYC vs. siControl (design `~ CellLine + Condition`). |
| `03b_erk_deseq2.R` | `salmon_quant_erk/`, `data/erk_metadata.csv` | `results/ERK_DE_<tp>_vs_vehicle.csv` | DESeq2: each ERKi timepoint vs. vehicle (design `~ CellLine + Timepoint`). |
| `04_plot_gene_comparison.R` | both DESeq2 results | `results/gene_comparison_barplot.pdf` | log2FC bar plot for the 6-gene panel across siMYC + ERKi timepoints. |

