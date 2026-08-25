import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from _common import apply_gene_aliases

FIGDIR = 'results/figures'

adata = sc.read_h5ad('processed_annotated.h5ad')
apply_gene_aliases(adata)  # TSTA3 -> GFUS
adata_mal = adata[adata.obs['Cell_type'] == 'Malignant ductal cells'].copy()
GENE_PANEL = ['MYC', 'HAS2', 'GPAA1', 'ST3GAL1', 'EXT1', 'B3GNT3', 'GFUS']
present = [g for g in GENE_PANEL if g in adata_mal.var_names]

raw = adata_mal[:, present].layers['counts']
if not isinstance(raw, np.ndarray):
    raw = raw.toarray()
count_df = pd.DataFrame(raw, columns=present, index=adata_mal.obs_names)

# Top 2 GTs by Jaccard with MYC, pinned to the values computed against the
# published dataset (GFUS: 0.275, GPAA1: 0.254 — both ahead of B3GNT3: 0.228,
# ST3GAL1: 0.224, EXT1: 0.185, HAS2: 0.028) so the figure doesn't silently
# change genes if the input data changes. Set PIN_GT1_GT2 = False to
# re-derive them from the current data instead.
PIN_GT1_GT2 = True
GT1, GT2 = 'GFUS', 'GPAA1'

if not PIN_GT1_GT2:
    myc_on = count_df['MYC'] > 0
    jaccard_scores = {}
    for g in [x for x in present if x != 'MYC']:
        gt_on = count_df[g] > 0
        both = (myc_on & gt_on).sum()
        either = (myc_on | gt_on).sum()
        jaccard_scores[g] = both / either if either > 0 else 0
    GT1, GT2 = pd.Series(jaccard_scores).sort_values(ascending=False).head(2).index.tolist()

print(f'Top 2 GTs: {GT1}, {GT2}')

categories = ['none', 'MYC only', f'{GT1} only', f'{GT2} only',
              f'MYC+{GT1}', f'MYC+{GT2}', f'{GT1}+{GT2}', 'all three']

def assign_bin(myc, gt1, gt2):
    if myc and gt1 and gt2: return 'all three'
    elif myc and gt1: return f'MYC+{GT1}'
    elif myc and gt2: return f'MYC+{GT2}'
    elif gt1 and gt2: return f'{GT1}+{GT2}'
    elif myc: return 'MYC only'
    elif gt1: return f'{GT1} only'
    elif gt2: return f'{GT2} only'
    else: return 'none'

records = []
for t in range(6):
    myc_b, gt1_b, gt2_b = count_df['MYC'] > t, count_df[GT1] > t, count_df[GT2] > t
    bins = [assign_bin(m, g1, g2) for m, g1, g2 in zip(myc_b, gt1_b, gt2_b)]
    bc = pd.Series(bins).value_counts(normalize=True)
    for cat in categories:
        records.append({'threshold': t, 'category': cat, 'fraction': bc.get(cat, 0)})

sweep_df = pd.DataFrame(records)
print(sweep_df.pivot(index='threshold', columns='category', values='fraction'))

fig, ax = plt.subplots(figsize=(8, 6))
for cat in categories:
    sub = sweep_df[sweep_df['category'] == cat]
    ax.plot(sub['threshold'], sub['fraction'], marker='o', label=cat)
ax.axvline(x=1, linestyle=':', color='grey')
ax.set_xlabel('Detection threshold (count > x)')
ax.set_ylabel('Fraction of malignant cells')
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.savefig(f'{FIGDIR}/glyco_panel_sensitivity_sweep.pdf', bbox_inches='tight')
