import os
import scanpy as sc
from _common import move_category_first

IN_H5AD = 'adata_relabeled.h5ad'  # produced by 01b_quick_relabel.py
TARGET = 'Malignant ductal cells'

# --- Load the shared, already-normalized/relabeled object ---
adata = sc.read_h5ad(IN_H5AD)

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

# --- Move malignant group to front of the groupby axis ---
move_category_first(adata_t, 'Cell_type', TARGET)

# --- Rank genes per cell type (tumor cells only) to get pts specificity ---
# This is the slow step (wilcoxon across every gene) -- cache the result below
# so plotting scripts never have to rerun it.
sc.tl.rank_genes_groups(adata_t, groupby='Cell_type', method='wilcoxon', pts=True)
pts = adata_t.uns['rank_genes_groups']['pts']

# --- Cache for downstream plotting scripts ---
os.makedirs('cache', exist_ok=True)
adata_t.write('cache/adata_t_processed.h5ad')
pts.to_pickle('cache/pts.pkl')

print('Done. Cached adata_t and pts to cache/')
print('Rerun this script only when the source data or the ranking logic changes.')
