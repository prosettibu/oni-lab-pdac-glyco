import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from _common import move_category_first, specificity_order

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR
sc.settings.verbosity = 2

TARGET = 'Malignant ductal cells'

# --- Load the shared, already-normalized/relabeled object (see 01b_quick_relabel.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

# =========================================================================
# Core glycosyltransferase panel (the only genes used)
# EXT1, HAS2, ST3GAL1, and MYC are excluded entirely — none of these were
# in the original list, so they're left out regardless of DE results.
# =========================================================================
genes = ['NDRG1', 'MAL2', 'ENPP2', 'PTK2', 'CDH17', 'PLEC', 'EIF3E', 'SDCBP', 'SULF1',
         'SDC2', 'LY6D', 'MAPK15', 'PSCA', 'SCRIB', 'SLURP1', 'SLURP2', 'THEM6',
         'VPS28', 'AGO2', 'CCN3', 'COL14A1', 'RAB2A', 'PXDNL']

# --- Tumor only, for the actual panel plots ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()
move_category_first(adata_t, 'Cell_type', TARGET)

specificity = specificity_order(adata_t, groupby='Cell_type', target=TARGET, genes=genes)
present_ordered = specificity.index.tolist()
print(specificity)

# --- Violin plots: one gene per PDF page, tumor cells only, uniqueness order ---
with PdfPages(f'{FIGDIR}/violin_glyco_panel_tumor_by_celltype.pdf') as pdf:
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
              colorbar_title='Relative mean expression\n(scaled per gene, 0–1)',
              size_title='Fraction of cells\nexpressing (%)',
              save='_glyco_panel_tumor_by_celltype.pdf')

print(f'Done. Figures in {FIGDIR}/')
