library(ggplot2)
library(dplyr)

genes <- c("EXT1", "ST3GAL1", "GPAA1", "HAS2", "MYC", "GFUS")

# Load siMYC results (see scripts/03_simyc_deseq2.R)
simyc <- read.csv("results/siMYC_DE_results.csv")
simyc <- simyc %>% filter(gene_name %in% genes) %>%
  mutate(Condition = "siMYC") %>%
  select(gene_name, log2FoldChange, padj, Condition)

# Load ERK time course results (see scripts/03b_erk_deseq2.R)
erk_list <- lapply(c("1h", "4h", "12h", "24h"), function(tp) {
  df <- read.csv(paste0("results/ERK_DE_", tp, "_vs_vehicle.csv"))
  df %>% filter(gene_name %in% genes) %>%
    mutate(Condition = paste0("ERKi_", tp)) %>%
    select(gene_name, log2FoldChange, padj, Condition)
})
erk <- bind_rows(erk_list)

# Combine
plot_data <- bind_rows(simyc, erk)
plot_data$Condition <- factor(plot_data$Condition,
                               levels = c("siMYC", "ERKi_1h", "ERKi_4h", "ERKi_12h", "ERKi_24h"))
plot_data$Significant <- ifelse(!is.na(plot_data$padj) & plot_data$padj < 0.05, "padj < 0.05", "n.s.")

p <- ggplot(plot_data, aes(x = gene_name, y = log2FoldChange, fill = Significant)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~Condition, nrow = 1) +
  scale_fill_manual(values = c("padj < 0.05" = "firebrick", "n.s." = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Gene", y = "log2 Fold Change", title = "siMYC vs ERKi time course: 5 genes of interest")

dir.create("results", showWarnings = FALSE)
ggsave("results/gene_comparison_barplot.pdf", p, width = 12, height = 5)
cat("Saved plot to results/gene_comparison_barplot.pdf\n")
