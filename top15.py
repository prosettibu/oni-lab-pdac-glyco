import scanpy as sc
import pandas as pd

sc.settings.figdir = './figures'
sc.settings.verbosity = 2

# --- Load cached, pre-normalized, relabeled object (see prep.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

target = 'Malignant ductal cells'

# --- Top 15 significant genes (cBioPortal Altered vs Unaltered) + B3GNT3 ---
genes = ['FOXL1', 'NUDT16L2P', 'B3GNT3', 'SLC35F2', 'SDC4', 'PLEK2', 'SLC26A11',
         'MAL2', 'STEAP3', 'IL1RN', 'ESYT3', 'LRRC1', 'COL17A1', 'TRIM29', 'ZFP3']

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

# --- Compute pts (fraction of cells expressing) per cell type ---
sc.tl.rank_genes_groups(adata_t, groupby='Cell_type', method='wilcoxon', pts=True)

# --- Order genes left-to-right by fraction of malignant ductal cells expressing ---
pts = adata_t.uns['rank_genes_groups']['pts']
malignant_pts = pts[target].loc[present].sort_values(ascending=False)
present_ordered = malignant_pts.index.tolist()

print(malignant_pts)

# --- Dot plot ---
sc.pl.dotplot(adata_t, present_ordered, groupby='Cell_type',
              standard_scale='var',
              colorbar_title='Relative mean expression\n(scaled per gene, 0\u20131)',
              size_title='Fraction of cells\nexpressing (%)',
              save='_cbioportal_top15_tumor_by_celltype.pdf')

print("Done. Figures in ./figures/")