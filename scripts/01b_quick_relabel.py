import scanpy as sc

# --- Load, normalize, relabel ONCE. Downstream scripts load the cached ---
# --- output of this instead of repeating normalize/log1p/relabel each time. ---

adata = sc.read_h5ad('StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
adata.var_names_make_unique()
adata.var.index.name = None  # avoid index/column name collision on write

adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(
    {'Ductal cell type 2': 'Malignant ductal cells', 'Fibroblast cell': 'Fibroblast', 'Macrophage cell': 'Macrophage'}
)

adata.write('adata_relabeled.h5ad')
print("Saved adata_relabeled.h5ad — downstream scripts (b3.py, glyco.py, myc.py, check.py) load this directly.")