import os
import scanpy as sc

sc.settings.figdir = './newfigures'
sc.settings.verbosity = 2

os.makedirs('./cache', exist_ok=True)

# --- Load and normalize ---
adata = sc.read_h5ad('StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
adata.var_names_make_unique()
adata.var.index.name = None

# --- Relabel malignant population ---
adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(
    {'Ductal cell type 2': 'Malignant ductal cells', 'Macrophage cell': 'Macrophage', 'Fibroblast cell': 'Fibroblast'}
)

target = 'Malignant ductal cells'

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

# --- Move malignant group to front of the groupby axis ---
adata_t.obs['Cell_type'] = adata_t.obs['Cell_type'].cat.reorder_categories(
    [target] + [c for c in adata_t.obs['Cell_type'].cat.categories if c != target]
)

# --- Rank genes per cell type (tumor cells only) to get pts specificity ---
# This is the slow step (wilcoxon across every gene) -- cache the result below
# so plotting scripts never have to rerun it.
sc.tl.rank_genes_groups(adata_t, groupby='Cell_type', method='wilcoxon', pts=True)
pts = adata_t.uns['rank_genes_groups']['pts']

# --- Cache for downstream plotting scripts ---
adata_t.write('./cache/adata_t_processed.h5ad')
pts.to_pickle('./cache/pts.pkl')

print('Done. Cached adata_t and pts to ./cache/')
print('Rerun this script only when the source data or the ranking logic changes.')