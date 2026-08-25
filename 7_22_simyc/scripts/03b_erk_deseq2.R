library(tximport)
library(DESeq2)

# Reference tx2gene mapping (GENCODE v46 / GRCh38) — not included in this
# repo. Point TX2GENE_PATH at your own copy, or set the env var.
tx2gene_path <- Sys.getenv("TX2GENE_PATH", "../reference/human/tx2gene.tsv")
tx2gene <- read.table(tx2gene_path,
                       header = FALSE, sep = "\t",
                       col.names = c("transcript_id", "gene_id", "gene_name"))

# Load sample metadata
meta <- read.csv("data/erk_metadata.csv")
meta$Timepoint <- factor(meta$Timepoint, levels = c("vehicle", "1h", "4h", "12h", "24h"))
meta$CellLine <- factor(meta$CellLine)

# Build paths to salmon quant.sf files (see scripts/02b_salmon_quant_erk.sh)
quant_dir <- "salmon_quant_erk"
files <- file.path(quant_dir, meta$Run, "quant.sf")
names(files) <- meta$Run

stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                 ignoreAfterBar = TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ CellLine + Timepoint)
dds <- DESeq(dds)

gene_map <- unique(tx2gene[, c("gene_id", "gene_name")])

dir.create("results", showWarnings = FALSE)

# Extract results for each timepoint vs vehicle
for (tp in c("1h", "4h", "12h", "24h")) {
  res <- results(dds, contrast = c("Timepoint", tp, "vehicle"))
  res <- as.data.frame(res)
  res$gene_id <- rownames(res)
  res <- merge(res, gene_map, by = "gene_id", all.x = TRUE)
  outfile <- paste0("results/ERK_DE_", tp, "_vs_vehicle.csv")
  write.csv(res, outfile, row.names = FALSE)
  cat("Wrote", outfile, "\n")
}

cat("Done.\n")
