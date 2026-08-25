import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import scanpy as sc
from _common import apply_gene_aliases

FIGDIR = 'results/figures'

# --- Load the QC/clustered object (see 01_qc_clustering_annotation.py) —
# --- needed here for its raw 'counts' layer, which the lighter
# --- adata_relabeled.h5ad (01b) doesn't carry. ---
adata = sc.read_h5ad('processed_annotated.h5ad')
apply_gene_aliases(adata)  # TSTA3 -> GFUS
adata_mal = adata[adata.obs['Cell_type'] == 'Malignant ductal cells'].copy()

GENE_PANEL_NO_B3GNT3 = ['MYC', 'HAS2', 'GPAA1', 'ST3GAL1', 'EXT1', 'GFUS']
present_no_b3gnt3 = [g for g in GENE_PANEL_NO_B3GNT3 if g in adata_mal.var_names]

raw = adata_mal[:, present_no_b3gnt3].layers['counts']
if not isinstance(raw, np.ndarray):
    raw = raw.toarray()
binary_df = pd.DataFrame(raw > 0, columns=present_no_b3gnt3, index=adata_mal.obs_names)

bool_df = binary_df[present_no_b3gnt3].astype(bool)
bool_df = bool_df[bool_df.any(axis=1)]  # drop cells with nothing on

data = from_indicators(present_no_b3gnt3, bool_df)
upset = UpSet(data, subset_size='count', sort_by='cardinality',
              max_subset_rank=25, show_counts=False)
fig = plt.figure(figsize=(12, 6))
upset.plot(fig=fig)
fig.suptitle('MYC / glycotransferase co-detection (no B3GNT3) — malignant ductal cells')
plt.savefig(f'{FIGDIR}/glyco_myc_upset_no_b3gnt3.pdf', bbox_inches='tight')
plt.close('all')
