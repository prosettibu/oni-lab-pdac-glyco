import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sc.settings.figdir = './figures'

# --- Load and normalize ---
adata = sc.read_h5ad('StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
adata.var_names_make_unique()

# --- Relabel malignant population ---
adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(
    {'Ductal cell type 2': 'Malignant ductal cells'}
)

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

genes = ['GALNT4','GALNT6','GALNT7','MGAT3','MGAT4A','MGAT4B','MGAT4C','MGAT5',
         'FUT1','FUT2','FUT4','FUT8','ST3GAL6','ST6GALNAC2','ST6GALNAC4',
         'EXT1','HAS2','ST3GAL1',
         'B3GNT2','B3GNT3','B4GALT4','B4GALT5','B4GALNT2','B4GALNT3','MYC',
         'B3GALT4','B3GALT5','ABO','GGTA1','GCNT1','GCNT3','GCNT4','RFNG',
         'EXTL1','CHPF','UGCG','POMT1','COLGALT2','HAS3','ALG3','ALG13','PIGZ']

present = [g for g in genes if g in adata_t.var_names]
missing = [g for g in genes if g not in adata_t.var_names]
print('Found:', len(present), present)
if missing:
    print('Missing from dataset:', missing)

# --- Rank genes per cell type (tumor cells only) to get pts specificity ---
sc.tl.rank_genes_groups(adata_t, groupby='Cell_type', method='wilcoxon', pts=True)

target = 'Malignant ductal cells'
pts = adata_t.uns['rank_genes_groups']['pts']  # genes x groups, fraction expressing

specificity = pts[target] - pts.drop(columns=target).max(axis=1)
specificity = specificity.loc[present].sort_values(ascending=False)
present_ordered = specificity.index.tolist()

print(specificity)  # sanity check

# --- Move malignant group to front of the groupby axis ---
adata_t.obs['Cell_type'] = adata_t.obs['Cell_type'].cat.reorder_categories(
    [target] + [c for c in adata_t.obs['Cell_type'].cat.categories if c != target]
)

# --- Violin plots: one gene per PDF page, tumor cells only, by cell type, uniqueness order ---
with PdfPages('./figures/violin_NvsM_upregM_tumor_by_celltype.pdf') as pdf:
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

print("Done. PDFs saved in ./figures/")