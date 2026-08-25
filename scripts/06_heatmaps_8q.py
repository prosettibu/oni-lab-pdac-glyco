import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FIGDIR = 'results/figures'
sc.settings.figdir = FIGDIR
os.makedirs(FIGDIR, exist_ok=True)

# --- Load cached data (produced by 02_compute_specificity_stats.py) ---
adata_t = sc.read_h5ad('cache/adata_t_processed.h5ad')
pts = pd.read_pickle('cache/pts.pkl')
target = 'Malignant ductal cells'

# --- Color scheme: navy blue -> neon green ---
navy_neon = LinearSegmentedColormap.from_list('navy_neon', ['#00004D', '#39FF14'])

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
    plt.savefig(f'{FIGDIR}/{save_name}', bbox_inches='tight')
    plt.close('all')


# --- Get all genes on 8q via Ensembl BioMart, cached locally after the first
# --- run so reruns don't depend on network access or on Ensembl's annotation
# --- content staying identical over time. ---
GENES_8Q_CACHE = 'cache/genes_8q_biomart.csv'

if os.path.exists(GENES_8Q_CACHE):
    genes_8q = pd.read_csv(GENES_8Q_CACHE)['gene'].tolist()
    print(f'Loaded {len(genes_8q)} chr8q genes from {GENES_8Q_CACHE}')
else:
    biomart = sc.queries.biomart_annotations(
        'hsapiens',
        ['external_gene_name', 'chromosome_name', 'band']
    ).query("chromosome_name == '8' and band.str.startswith('q')", engine='python')

    genes_8q = biomart['external_gene_name'].dropna().unique().tolist()
    os.makedirs('cache', exist_ok=True)
    pd.DataFrame({'gene': genes_8q}).to_csv(GENES_8Q_CACHE, index=False)
    print(f'Queried BioMart: {len(genes_8q)} chr8q genes, cached to {GENES_8Q_CACHE}')

plot_percell_heatmap(genes_8q, 'heatmap_8q_percell_vertical.pdf')

print(f'Done. Figures in {FIGDIR}/')
