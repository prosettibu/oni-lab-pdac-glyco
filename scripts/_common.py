"""Shared helpers for the 7_6_scrnaseq pipeline scripts.

Kept deliberately small: a gene-symbol alias table, the canonical
cell-type relabeling used by every script, and the two computations
(reorder-target-first, specificity-vs-other-groups) that were being
copy-pasted across the panel-plotting scripts.
"""
import scanpy as sc

# TSTA3 and GFUS are the same gene (alternate symbols). The source data
# uses TSTA3; the papers/panels this pipeline reproduces use GFUS.
GENE_ALIASES = {'TSTA3': 'GFUS'}

CELL_TYPE_RENAME = {
    'Ductal cell type 2': 'Malignant ductal cells',
    'Macrophage cell': 'Macrophage',
    'Fibroblast cell': 'Fibroblast',
}


def apply_gene_aliases(adata, aliases=GENE_ALIASES):
    for old, new in aliases.items():
        adata.var_names = adata.var_names.where(adata.var_names != old, new)
    return adata


def relabel_cell_types(adata, rename=CELL_TYPE_RENAME):
    adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(rename)
    return adata


def move_category_first(adata, obs_col, target):
    """Reorder `obs_col`'s categories so `target` plots first (e.g. left column of a dotplot)."""
    cats = adata.obs[obs_col].cat.categories
    adata.obs[obs_col] = adata.obs[obs_col].cat.reorder_categories(
        [target] + [c for c in cats if c != target]
    )
    return adata


def specificity_order(adata, groupby, target, genes):
    """Rank `genes` (present in var_names) by how specific their detection is to
    `target` vs. the next-highest-detecting other group in `groupby`.

    Runs sc.tl.rank_genes_groups(method='wilcoxon', pts=True) as a side effect
    and returns a Series of genes present in the data, sorted most- to least-specific.
    """
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    print('Found:', len(present), present)
    if missing:
        print('Missing from dataset:', missing)

    sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon', pts=True)
    pts = adata.uns['rank_genes_groups']['pts']
    specificity = pts[target].loc[present] - pts.drop(columns=target).loc[present].max(axis=1)
    return specificity.sort_values(ascending=False)
