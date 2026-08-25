"""WIP: immune/myeloid marker dotplot. Not part of the numbered pipeline —
kept here pending review before promotion into the main sequence.
"""
import scanpy as sc
from _common import move_category_first

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR
sc.settings.verbosity = 2

TARGET = 'Malignant ductal cells'

# --- Load the shared, already-normalized/relabeled object (see 01b_quick_relabel.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()
move_category_first(adata_t, 'Cell_type', TARGET)

# =========================================================================
# Immune/myeloid marker panel
# =========================================================================
genes = ['IL4', 'IL4R', 'IL2RG', 'CLEC10A', 'TREM1', 'CCL2', 'CDKN1A',
         'CCL7', 'CCR2', 'CCR1', 'CCR3', 'CCR5', 'GPNMB']
present = [g for g in genes if g in adata_t.var_names]
missing = [g for g in genes if g not in adata_t.var_names]
print('Found:', len(present), present)
if missing:
    print('Missing from dataset:', missing)

sc.pl.dotplot(adata_t, present, groupby='Cell_type',
              standard_scale='var',
              colorbar_title='Relative mean expression\n(scaled per gene, 0–1)',
              size_title='Fraction of cells\nexpressing (%)',
              save='_charlotte_tumor_by_celltype.pdf')
