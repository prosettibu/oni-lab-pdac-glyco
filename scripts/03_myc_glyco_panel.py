import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

navy_neon = LinearSegmentedColormap.from_list('navy_neon', ['#000080', '#39FF14'])

sc.settings.figdir = './newfigures'
sc.settings.verbosity = 2

# --- Load and normalize ---
adata = sc.read_h5ad('StdWf1_PRJCA001063_CRC_besca2.raw.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.var_names = adata.var['SYMBOL'].astype(str)
adata.var_names = adata.var_names.where(adata.var_names != 'TSTA3', 'GFUS')
adata.var_names_make_unique()

# --- Relabel malignant population ---
adata.obs['Cell_type'] = adata.obs['Cell_type'].cat.rename_categories(
    {'Ductal cell type 2': 'Malignant ductal cells', 'Macrophage cell': 'Macrophage', 'Fibroblast cell': 'Fibroblast'}
)

target = 'Malignant ductal cells'

# --- MYC-paper gene set ---
genes = ['EXT1', 'HAS2', 'ST3GAL1', 'MYC', 'B3GNT3', 'GPAA1', 'GFUS']

# --- Tumor only ---
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

# =========================================================================
# DOT PLOT
# =========================================================================
sc.pl.dotplot(adata_t, present_ordered, groupby='Cell_type',
              standard_scale='var',
              colorbar_title='Relative mean expression\n(scaled per gene, 0\u20131)',
              size_title='Fraction of cells\nexpressing (%)',
              save='_myc_panel_tumor_by_celltype.pdf')

# =========================================================================
# VIOLIN PLOTS (one gene per PDF page, uniqueness order)
# =========================================================================
with PdfPages('./figures/violin_myc_panel_tumor_by_celltype.pdf') as pdf:
    for gene in present_ordered:
        with plt.rc_context({'figure.figsize': (8, 8)}):
            sc.pl.violin(adata_t, gene, groupby='Cell_type', rotation=90,
                         stripplot=False, show=False)
            fig = plt.gcf()
            ax = fig.axes[0]
            ax.set_ylabel('Log-normalized expression')
            ax.set_title(gene)
            fig.subplots_adjust(bottom=0.45)
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

# =========================================================================
# HEATMAP OPTION A: per-cell expression (sc.pl.heatmap), vertical orientation
# =========================================================================
def plot_percell_heatmap(gene_list, save_name, cmap=navy_neon,
                          cbar_title='Log-normalized\nexpression'):
    present = [g for g in gene_list if g in adata_t.var_names]
    missing = [g for g in gene_list if g not in adata_t.var_names]
    if missing:
        print(f'Missing from dataset ({save_name}):', missing)

    spec = pts[target].loc[present] - pts.drop(columns=target).loc[present].max(axis=1)
    ordered = spec.sort_values(ascending=False).index.tolist()

    axes_dict = sc.pl.heatmap(
        adata_t, ordered, groupby='Cell_type',
        swap_axes=False, cmap=cmap, show=False, save=False
    )
    cbar_ax = axes_dict.get('color_legend_ax') or plt.gcf().axes[-1]
    cbar_ax.set_title(cbar_title, fontsize=10, pad=8)
    plt.savefig(f'./figures/{save_name}', bbox_inches='tight')
    plt.close('all')

genes_panel1 = ['EXT1', 'HAS2', 'MYC', 'ST3GAL1', 'GPAA1','GFUS']
genes_panel2 = ['NDRG1', 'MAL2', 'ENPP2', 'PTK2', 'CDH17', 'PLEC', 'EIF3E',
                'SDCBP', 'SULF1', 'SDC2', 'LY6D', 'MAPK15', 'PSCA', 'SCRIB',
                'SLURP1', 'THEM6', 'VPS28', 'AGO2', 'CCN3', 'COL14A1',
                'RAB2A', 'PXDNL']

plot_percell_heatmap(genes_panel1, 'heatmap_myc_panel_percell_vertical.pdf')   # replaces existing file
plot_percell_heatmap(genes_panel2, 'heatmap_panel2_percell_vertical.pdf')      # new panel

# =========================================================================
# HEATMAP OPTION B: aggregated mean expression per cell type, min-max scaled
# 0-1 per gene, cividis colormap, vertical layout, larger colorbar.
# =========================================================================
expr = adata_t[:, present_ordered].X
expr = np.asarray(expr.todense()) if hasattr(expr, 'todense') else np.asarray(expr)

df = pd.DataFrame(expr, columns=present_ordered)
df['Cell_type'] = adata_t.obs['Cell_type'].values
mean_expr = df.groupby('Cell_type')[present_ordered].mean()

cat_order = adata_t.obs['Cell_type'].cat.categories
mean_expr = mean_expr.reindex(cat_order)

# min-max scale each gene (column) independently to 0-1
scaled = (mean_expr - mean_expr.min()) / (mean_expr.max() - mean_expr.min())

fig, ax = plt.subplots(figsize=(1.2 * len(present_ordered) + 1.5, 6))
im = ax.imshow(scaled.values, cmap='cividis', aspect='auto', vmin=0, vmax=1)

ax.set_yticks(range(len(scaled.index)))
ax.set_yticklabels(scaled.index)
ax.set_xticks(range(len(present_ordered)))
ax.set_xticklabels(present_ordered, rotation=90)

ax.set_yticks(np.arange(-0.5, len(scaled.index), 1), minor=True)
ax.set_xticks(np.arange(-0.5, len(present_ordered), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.5)
ax.tick_params(which='minor', length=0)

cbar = fig.colorbar(im, ax=ax, fraction=0.1, pad=0.05)
cbar.ax.tick_params(labelsize=14)
cbar.outline.set_linewidth(1.2)
cbar.set_label('Relative mean expression\n(scaled per gene, 0\u20131)')

plt.tight_layout()
plt.savefig('./figures/heatmap_myc_panel_aggregated_scaled.pdf')
plt.close()

print("Done. Figures in ./figures/:")
print(" - dotplot_myc_panel_tumor_by_celltype.pdf")
print(" - violin_myc_panel_tumor_by_celltype.pdf")
print(" - heatmap_myc_panel_percell_vertical.pdf   (Option A: per-cell strips)")
print(" - heatmap_myc_panel_aggregated_scaled.pdf  (Option B: aggregated mean, 0-1 scaled)")