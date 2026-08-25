library(tximport)
library(DESeq2)

# Reference tx2gene mapping (GENCODE v46 / GRCh38) — not included in this
# repo. Point TX2GENE_PATH at your own copy, or set the env var.
tx2gene_path <- Sys.getenv("TX2GENE_PATH", "../reference/human/tx2gene.tsv")
tx2gene <- read.table(tx2gene_path,
                       header = FALSE, sep = "\t",
                       col.names = c("transcript_id", "gene_id", "gene_name"))

# Load sample metadata
meta <- read.csv("data/simyc_metadata.csv")
meta$Condition <- factor(meta$Condition, levels = c("siControl", "siMYC"))

# Build paths to salmon quant.sf files (see scripts/02_salmon_quant_simyc.sh)
quant_dir <- "salmon_quant_simyc"
files <- file.path(quant_dir, meta$Run, "quant.sf")
names(files) <- meta$Run

# Confirm all files exist before proceeding
stopifnot(all(file.exists(files)))

# Import and collapse transcript-level counts to gene-level
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                 ignoreAfterBar = TRUE)

# Build DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ CellLine + Condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "siMYC", "siControl"))
res <- as.data.frame(res)
res$gene_id <- rownames(res)

# Attach gene symbols
gene_map <- unique(tx2gene[, c("gene_id", "gene_name")])
res <- merge(res, gene_map, by = "gene_id", all.x = TRUE)

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/siMYC_DE_results.csv", row.names = FALSE)

cat("Done. Results written to results/siMYC_DE_results.csv\n")
