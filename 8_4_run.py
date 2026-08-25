import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("processed_annotated.h5ad")
blocks = pd.read_csv("../organoid_bulk/results/hM1_8q_blocks.csv")

genes_present = blocks[blocks["symbol"].isin(adata.var_names)].copy()
missing = set(blocks["symbol"]) - set(genes_present["symbol"])
print(f"{len(missing)} 8q genes not found in scRNA-seq data")

adata.obs["CONDITION"] = adata.obs["CONDITION"].astype("category")
adata.obs["CONDITION"] = adata.obs["CONDITION"].cat.reorder_categories(
    ["hM1A", "hM1E", "hM1F"]  # ADJUST to match your actual values from the print() check
)

for block_name, block_genes in genes_present.groupby("block"):
    gene_list = sorted(block_genes["symbol"].tolist())
    if len(gene_list) == 0:
        continue

    sc.pl.dotplot(
        adata,
        var_names=gene_list,
        groupby="CONDITION",
        standard_scale="var",
        title=f"{block_name} (n={len(gene_list)} genes) - A to F trajectory",
        save=f"_8q_{block_name}_trajectory.pdf"
    )



