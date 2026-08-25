import scanpy as sc
from _common import apply_gene_aliases, relabel_cell_types

RAW_H5AD = 'StdWf1_PRJCA001063_CRC_besca2.raw.h5ad'
OUT_H5AD = 'adata_relabeled.h5ad'

# --- Load, normalize, relabel ONCE. Downstream scripts load OUT_H5AD instead ---
# --- of repeating normalize/log1p/relabel each time: 02_compute_specificity_stats.py, ---
# --- 03_myc_glyco_panel.py, 03b_glyco_pathway_panel.py, 03c_nglycan_oglycan_panel.py, ---
# --- 03d_b3gnt3_focused.py, 03e_cbioportal_top15.py, 99_wip_immune_panel_dotplot.py ---

adata = sc.read_h5ad(RAW_H5AD)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
apply_gene_aliases(adata)  # TSTA3 -> GFUS
adata.var_names_make_unique()
adata.var.index.name = None  # avoid index/column name collision on write

relabel_cell_types(adata)

adata.write(OUT_H5AD)
print(f'Saved {OUT_H5AD}')
