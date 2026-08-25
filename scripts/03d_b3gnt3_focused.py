import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _common import move_category_first

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR

# --- Load the shared, already-normalized/relabeled object (see 01b_quick_relabel.py) ---
adata = sc.read_h5ad('adata_relabeled.h5ad')

# --- Tumor only ---
adata_t = adata[adata.obs['CONDITION'] == 'T'].copy()

gene = 'B3GNT3'
assert gene in adata_t.var_names, f"{gene} not found in dataset"

# --- Move malignant group to front of the groupby axis ---
target = 'Malignant ductal cells'
move_category_first(adata_t, 'Cell_type', target)

# --- Dot plot ---
sc.pl.dotplot(adata_t, [gene], groupby='Cell_type',
              save='_B3GNT3_tumor_by_celltype.pdf')

# --- Violin plot (taller figure + extra bottom margin so labels don't clip) ---
with plt.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.violin(adata_t, gene, groupby='Cell_type', rotation=90,
                 stripplot=False, show=False)
    fig = plt.gcf()
    ax = fig.axes[0]
    ax.set_ylabel('Log-normalized expression')
    ax.set_title(gene)
    fig.subplots_adjust(bottom=0.45)
    fig.savefig(f'{FIGDIR}/violin_B3GNT3_tumor_by_celltype.pdf', bbox_inches='tight')
    plt.close(fig)

# =========================================================================
# HEATMAP OPTION A: per-cell expression (sc.pl.heatmap), vertical orientation
# Cell types run down the y-axis, single gene column, one strip per cell.
# =========================================================================
sc.pl.heatmap(
    adata_t, [gene], groupby='Cell_type',
    swap_axes=False,          # cell types on y-axis (vertical), gene on x-axis
    cmap='cividis',
    save='_B3GNT3_percell_vertical.pdf'
)
plt.close('all')

# =========================================================================
# HEATMAP OPTION B: aggregated mean expression per cell type, min-max scaled
# 0-1, cividis colormap, vertical layout, larger colorbar. Matches the
# earlier custom style built from the reference image.
# =========================================================================
expr = adata_t[:, gene].X
expr = np.asarray(expr.todense()).flatten() if hasattr(expr, 'todense') else np.asarray(expr).flatten()

df = pd.DataFrame({'Cell_type': adata_t.obs['Cell_type'].values, 'expr': expr})
mean_expr = df.groupby('Cell_type')['expr'].mean()

cat_order = adata_t.obs['Cell_type'].cat.categories
mean_expr = mean_expr.reindex(cat_order)

scaled = (mean_expr - mean_expr.min()) / (mean_expr.max() - mean_expr.min())

fig, ax = plt.subplots(figsize=(3, 6))
im = ax.imshow(scaled.values.reshape(-1, 1), cmap='cividis', aspect='auto', vmin=0, vmax=1)

ax.set_yticks(range(len(scaled)))
ax.set_yticklabels(scaled.index)
ax.set_xticks([0])
ax.set_xticklabels([gene], rotation=90)
ax.set_xlim(-0.5, 0.5)

ax.set_yticks(np.arange(-0.5, len(scaled), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.5)
ax.tick_params(which='minor', length=0)

# larger colorbar: bigger fraction (width) and larger tick/label fonts
cbar = fig.colorbar(im, ax=ax, fraction=0.25, pad=0.15)
cbar.ax.tick_params(labelsize=14)
cbar.outline.set_linewidth(1.2)

plt.tight_layout()
plt.savefig(f'{FIGDIR}/heatmap_B3GNT3_aggregated_scaled.pdf')
plt.close()

print(f"Done. Figures in {FIGDIR}/:")
print(" - dotplot_B3GNT3_tumor_by_celltype.pdf")
print(" - violin_B3GNT3_tumor_by_celltype.pdf")
print(" - heatmap_B3GNT3_percell_vertical.pdf   (Option A: per-cell strips)")
print(" - heatmap_B3GNT3_aggregated_scaled.pdf  (Option B: aggregated mean, 0-1 scaled)")