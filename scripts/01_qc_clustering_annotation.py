import scanpy as sc
import pandas as pd
from _common import apply_gene_aliases, relabel_cell_types, CELL_TYPE_RENAME

RAW_H5AD = 'StdWf1_PRJCA001063_CRC_besca2.raw.h5ad'
OUT_H5AD = 'processed_annotated.h5ad'
RANDOM_STATE = 0

sc.settings.figdir = 'results/figures'
sc.settings.verbosity = 2

# --- Load ---
adata = sc.read_h5ad(RAW_H5AD)

# --- Gene symbol aliases (TSTA3 -> GFUS), applied before any snapshot (layers,
# --- .raw) is taken so every downstream copy of var_names stays consistent ---
apply_gene_aliases(adata)

# --- Relabel malignant population (full map: see _common.CELL_TYPE_RENAME) ---
relabel_cell_types(adata)
print(adata.obs['Cell_type'].value_counts())

# --- QC overview ---
adata.var['mt'] = adata.var['SYMBOL'].str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# --- Normalize ---
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# --- Dimensionality reduction (seeded for reproducibility) ---
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50, random_state=RANDOM_STATE)
sc.pp.neighbors(adata, random_state=RANDOM_STATE)
sc.tl.umap(adata, random_state=RANDOM_STATE)
sc.pl.umap(adata, color='Cell_type', save='_celltype.pdf')

# --- Marker genes per cell type, with pts (fraction expressing) ---
sc.tl.rank_genes_groups(adata, groupby='Cell_type', method='wilcoxon', pts=True)

markers = {}
groups = adata.obs['Cell_type'].unique().tolist()
for group in groups:
    names = adata.uns['rank_genes_groups']['names'][group][:5]
    markers[group] = list(names)
    print(group, ':', list(names))

marker_genes = list(dict.fromkeys([g for genes in markers.values() for g in genes]))

# --- Order marker_genes by specificity to malignant ductal cells ---
target = 'Malignant ductal cells'
pts = adata.uns['rank_genes_groups']['pts']

specificity = pts[target] - pts.drop(columns=target).max(axis=1)
specificity = specificity.loc[marker_genes].sort_values(ascending=False)
marker_genes_ordered = specificity.index.tolist()

print(specificity)

# --- Move malignant group to front of the groupby axis ---
adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.reorder_categories(
    [target] + [c for c in adata.obs['Cell_type'].cat.categories if c != target]
)

# --- Dot plot: derived markers, ordered by malignant-cell specificity ---
sc.pl.dotplot(adata, marker_genes_ordered, groupby='Cell_type',
              save='_derived_markers.pdf', standard_scale='var')

# --- Dot plot: split by condition too ---
sc.pl.dotplot(adata, marker_genes_ordered, groupby=['Cell_type', 'CONDITION'],
              save='_derived_markers_by_condition.pdf', standard_scale='var')

adata.write(OUT_H5AD)
print(f'Done. Figures in results/figures/, processed object saved to {OUT_H5AD}.')
