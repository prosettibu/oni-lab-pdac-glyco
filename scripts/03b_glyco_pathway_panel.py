import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.figdir = './newfigures'
sc.settings.verbosity = 2

# --- Load cached, pre-normalized, relabeled object (see prep.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

target = 'Malignant ductal cells'

# =========================================================================
# STEP 1: Core glycosyltransferase panel (the only genes used)
# EXT1, HAS2, ST3GAL1, and MYC are excluded entirely — none of these were
# in the original list, so they're left out regardless of DE results.
# =========================================================================
genes = ['NDRG1','MAL2','ENPP2','PTK2','CDH17','PLEC','EIF3E','SDCBP','SULF1',
         'SDC2','LY6D','MAPK15','PSCA','SCRIB','SLURP1','SLURP2','THEM6',
         'VPS28','AGO2','CCN3','COL14A1','RAB2A','PXDNL']

# --- Tumor only, for the actual panel plots ---
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
present_ordered = specificity.index.tolist()

print(specificity)

# --- Violin plots: one gene per PDF page, tumor cells only, uniqueness order ---
with PdfPages('./figures/violin_glyco_panel_tumor_by_celltype.pdf') as pdf:
    for gene in present_ordered:
        with plt.rc_context({'figure.figsize': (8, 6)}):
            sc.pl.violin(adata_t, gene, groupby='Cell_type', rotation=90,
                         stripplot=False, show=False)
            fig = plt.gcf()
            ax = fig.axes[0]
            ax.set_ylabel('Log-normalized expression')
            fig.subplots_adjust(bottom=0.35)
            plt.title(gene)
            pdf.savefig(fig)
            plt.close(fig)

# --- Dot plot: all genes, tumor cells only, uniqueness order ---
sc.pl.dotplot(adata_t, present_ordered, groupby='Cell_type',
              standard_scale='var',
              colorbar_title='Relative mean expression\n(scaled per gene, 0\u20131)',
              size_title='Fraction of cells\nexpressing (%)',
              save='_glyco_panel_tumor_by_celltype.pdf')

print("Done. Figures in ./figures/")