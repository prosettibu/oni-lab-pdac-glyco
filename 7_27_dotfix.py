import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

navy_neon = LinearSegmentedColormap.from_list('navy_neon', ['#000080', '#39FF14'])

sc.settings.figdir = '/lab/solexa_oni/patrick/7_27_tasks'
sc.settings.verbosity = 2

# --- Load and normalize ---
adata = sc.read_h5ad('StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
adata.var_names = adata.var_names.where(adata.var_names != 'TSTA3', 'GFUS')
adata.var_names_make_unique()

# --- Relabel malignant population ---
adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(
    {'Ductal cell type 2': 'Malignant ductal cells', 'Macrophage cell': 'Macrophage', 'Fibroblast cell': 'Fibroblast'}
)

target = 'Malignant ductal cells'

# --- MYC-paper gene set ---
genes = ['EXT1', 'HAS2', 'ST3GAL1', 'GPAA1', 'MYC', 'GFUS']

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

present = [g for g in genes if g in adata_t.var_names]
missing = [g for g in genes if g not in adata_t.var_names]
print('Found:', len(present), present)
if missing:
    print('Missing from dataset:', missing)

# --- Move malignant group to front of the groupby axis ---
adata_t.obs['Cell_type'] = adata_t.obs['Cell_type'].cat.reorder_categories(
    [target] + [c for c in adata_t.obs['Cell_type'].cat.categories if c != target]
)

# --- Rank genes per cell type (tumor cells only) to get pts specificity ---
sc.tl.rank_genes_groups(adata_t, groupby='Cell_type', method='wilcoxon', pts=True)

pts = adata_t.uns['rank_genes_groups']['pts']
specificity = pts[target] - pts.drop(columns=target).max(axis=1)
specificity = specificity.loc[present].sort_values(ascending=False)

print(specificity)


# =========================================================================
# DOT PLOT 3: <fill in>
# =========================================================================
genes3 = ['IL4', 'IL4R','IL2RG', 'CLEC10A','TREM1','CCL2','CDKN1A', 'CCL7', 'CCR2', 'CCR1', 'CCR3', 'CCR5', 'GPNMB']  # <-- fill in gene names here
present3 = [g for g in genes3 if g in adata_t.var_names]
missing3 = [g for g in genes3 if g not in adata_t.var_names]
print('Found:', len(present3), present3)
if missing3:
    print('Missing from dataset:', missing3)

sc.pl.dotplot(adata_t, present3, groupby='Cell_type',
standard_scale='var',
colorbar_title='Relative mean expression\n(scaled per gene, 0\u20131)',
size_title='Fraction of cells\nexpressing (%)',
save='_charlotte_tumor_by_celltype.pdf')