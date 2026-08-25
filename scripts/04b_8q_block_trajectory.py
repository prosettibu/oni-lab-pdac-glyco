import scanpy as sc
import pandas as pd

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR

# 8q gene -> block assignments (from a separate organoid_bulk analysis, not
# part of this repo's scripts). Vendored copy: see data/external/README.md.
BLOCKS_CSV = 'data/external/hM1_8q_blocks.csv'

adata = sc.read_h5ad('processed_annotated.h5ad')
blocks = pd.read_csv(BLOCKS_CSV)

genes_present = blocks[blocks['symbol'].isin(adata.var_names)].copy()
missing = set(blocks['symbol']) - set(genes_present['symbol'])
print(f'{len(missing)} 8q genes not found in scRNA-seq data')

# NOTE: this dataset's CONDITION column is Tumor/Normal ('T'/'N'), not an
# hM1A/E/F organoid trajectory — the original hM1 trajectory ordering doesn't
# apply here. Grouping by CONDITION as-is instead.
adata.obs['CONDITION'] = adata.obs['CONDITION'].astype('category')

for block_name, block_genes in genes_present.groupby('block'):
    gene_list = sorted(block_genes['symbol'].tolist())
    if len(gene_list) == 0:
        continue

    sc.pl.dotplot(
        adata,
        var_names=gene_list,
        groupby='CONDITION',
        standard_scale='var',
        title=f'{block_name} (n={len(gene_list)} genes) — Tumor vs Normal',
        save=f'_8q_{block_name}_by_condition.pdf'
    )
