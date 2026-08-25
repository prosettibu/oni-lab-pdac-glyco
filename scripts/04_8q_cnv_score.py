import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

FIGDIR = 'results/figures'

# Per-cell chr8q CNV scores (from a separate inferCNV analysis, not part of
# this repo's scripts). Vendored copy: see data/external/README.md.
SCORES_CSV = 'data/external/infercnv_chr8q_scores.csv'

IMMUNE_LABELS = ['T cell', 'B cell', 'Macrophage cell']
MALIGNANT_LABEL = 'Malignant ductal cells'

scores = pd.read_csv(SCORES_CSV, index_col=0)

plot_df = scores[scores['Cell_type'].isin(IMMUNE_LABELS + [MALIGNANT_LABEL])].copy()
plot_df['Group'] = np.where(plot_df['Cell_type'] == MALIGNANT_LABEL, 'Malignant', 'Immune (reference)')

stat, p = mannwhitneyu(
    plot_df.loc[plot_df['Group'] == 'Malignant', 'chr8q_cnv_score'],
    plot_df.loc[plot_df['Group'] == 'Immune (reference)', 'chr8q_cnv_score']
)

fig, ax = plt.subplots(figsize=(5, 6))
sns.violinplot(data=plot_df, x='Group', y='chr8q_cnv_score', ax=ax,
                order=['Immune (reference)', 'Malignant'], inner='box')
ax.set_ylabel('8q (q-arm) CNV score')
ax.set_xlabel('')
ax.set_title(f'8q amplification signal — Mann-Whitney p = {p:.2e}')
plt.tight_layout()
plt.savefig(f'{FIGDIR}/chr8q_violin_immune_vs_malignant.pdf', bbox_inches='tight')
plt.close('all')

print(f'Mann-Whitney p = {p:.3e}')
print(plot_df.groupby('Group')['chr8q_cnv_score'].describe())
print(f'Saved to {FIGDIR}/chr8q_violin_immune_vs_malignant.pdf')

mal_scores = plot_df.loc[plot_df['Group'] == 'Malignant', 'chr8q_cnv_score']
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(mal_scores, bins=100, color='steelblue')
ax.set_xlabel('8q CNV score')
ax.set_ylabel('Number of malignant cells')
ax.set_title('Malignant cell 8q score distribution — looking for bimodality')
plt.savefig(f'{FIGDIR}/chr8q_malignant_histogram.pdf', bbox_inches='tight')
