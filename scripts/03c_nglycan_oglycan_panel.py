import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from _common import move_category_first, specificity_order

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR

TARGET = 'Malignant ductal cells'

# --- Load the shared, already-normalized/relabeled object (see 01b_quick_relabel.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

genes = ['GALNT4', 'GALNT6', 'GALNT7', 'MGAT3', 'MGAT4A', 'MGAT4B', 'MGAT4C', 'MGAT5',
         'FUT1', 'FUT2', 'FUT4', 'FUT8', 'ST3GAL6', 'ST6GALNAC2', 'ST6GALNAC4',
         'EXT1', 'HAS2', 'ST3GAL1',
         'B3GNT2', 'B3GNT3', 'B4GALT4', 'B4GALT5', 'B4GALNT2', 'B4GALNT3', 'MYC',
         'B3GALT4', 'B3GALT5', 'ABO', 'GGTA1', 'GCNT1', 'GCNT3', 'GCNT4', 'RFNG',
         'EXTL1', 'CHPF', 'UGCG', 'POMT1', 'COLGALT2', 'HAS3', 'ALG3', 'ALG13', 'PIGZ']

specificity = specificity_order(adata_t, groupby='Cell_type', target=TARGET, genes=genes)
present_ordered = specificity.index.tolist()
print(specificity)  # sanity check

# --- Move malignant group to front of the groupby axis ---
move_category_first(adata_t, 'Cell_type', TARGET)

# --- Violin plots: one gene per PDF page, tumor cells only, by cell type, uniqueness order ---
with PdfPages(f'{FIGDIR}/violin_NvsM_upregM_tumor_by_celltype.pdf') as pdf:
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

# --- Dot plot: all genes, tumor cells only, by cell type, uniqueness order ---
sc.pl.dotplot(adata_t, present_ordered, groupby='Cell_type',
              save='_NvsM_upregM_tumor_by_celltype.pdf')

print(f'Done. PDFs saved in {FIGDIR}/')
