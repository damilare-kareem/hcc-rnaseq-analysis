# 03_differential_expression.R
# Paired differential expression: HCC tumour vs matched non-tumour liver (DESeq2)

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
})

# Use several CPU cores to speed up the slow steps (set n_workers <- 1 if you get a worker error)
n_workers <- min(4, max(1, parallel::detectCores() - 1))
bp <- if (n_workers > 1) SnowParam(workers = n_workers) else SerialParam()

# ---------- 1. Fit the model (the slow step) ----------
# The fitted model is saved, so re-running this script skips the fit.
# Delete data/processed/dds_fitted.rds if you re-run Bit 2.
fitted_file <- "data/processed/dds_fitted.rds"
if (file.exists(fitted_file)) {
  dds <- readRDS(fitted_file)
  cat("Loaded the saved fitted model.\n")
} else {
  dds <- readRDS("data/processed/dds.rds")
  cat("Fitting DESeq2 on", n_workers, "core(s), started", format(Sys.time(), "%H:%M"), "\n")
  dds <- DESeq(dds, parallel = n_workers > 1, BPPARAM = bp)
  saveRDS(dds, fitted_file)
  cat("Model fitted at", format(Sys.time(), "%H:%M"), "\n")
}
coef_name <- "tissue_tumor_vs_non_tumor"

# Statistics (p-values, FDR) from the standard Wald test
res <- results(dds, name = coef_name, alpha = 0.05)

# Shrunken fold changes: more reliable for ranking and plotting, especially for low-count genes
cat("Shrinking fold changes, started", format(Sys.time(), "%H:%M"), "\n")
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm",
                        parallel = n_workers > 1, BPPARAM = bp)

# ---------- 2. Results table ----------
padj_cutoff <- 0.05
lfc_cutoff  <- 1      # log2 fold change of 1 = two-fold

de <- data.frame(
  gene_key        = rownames(res),
  symbol          = rowData(dds)$symbol,
  ensembl         = rowData(dds)$ensembl,
  baseMean        = res$baseMean,
  log2FC          = res_shrunk$log2FoldChange,
  log2FC_unshrunk = res$log2FoldChange,
  lfcSE           = res$lfcSE,
  stat            = res$stat,
  pvalue          = res$pvalue,
  padj            = res$padj
)
de$direction <- "Not significant"
de$direction[!is.na(de$padj) & de$padj < padj_cutoff & de$log2FC >=  lfc_cutoff] <- "Up in tumour"
de$direction[!is.na(de$padj) & de$padj < padj_cutoff & de$log2FC <= -lfc_cutoff] <- "Down in tumour"
de <- de[order(de$padj, -abs(de$log2FC)), ]

write.csv(de, "results/deseq2_all_genes.csv", row.names = FALSE)
write.csv(de[de$direction != "Not significant", ], "results/deseq2_significant_genes.csv", row.names = FALSE)
saveRDS(de, "data/processed/de_results.rds")

# ---------- 3. Volcano plot ----------
dir_cols <- c("Up in tumour" = "#D55E00", "Down in tumour" = "#0072B2", "Not significant" = "grey75")
plot_df  <- de[!is.na(de$padj), ]
plot_df$neglog10_padj <- -log10(pmax(plot_df$padj, 1e-300))   # FDR can underflow to 0; cap it

n_up   <- sum(de$direction == "Up in tumour")
n_down <- sum(de$direction == "Down in tumour")

top_up   <- head(plot_df[plot_df$direction == "Up in tumour", ], 10)
top_down <- head(plot_df[plot_df$direction == "Down in tumour", ], 10)
labels   <- rbind(top_up, top_down)

volcano <- ggplot(plot_df, aes(log2FC, neglog10_padj, colour = direction)) +
  geom_point(size = 0.9, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_text_repel(data = labels, aes(label = symbol), size = 3, colour = "black",
                  max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0) +
  scale_colour_manual(values = dir_cols,
                      labels = c("Up in tumour"   = paste0("Up in tumour (", n_up, ")"),
                                 "Down in tumour" = paste0("Down in tumour (", n_down, ")"),
                                 "Not significant" = "Not significant")) +
  labs(x = "log2 fold change (tumour vs matched non-tumour liver)",
       y = expression(-log[10]~"(adjusted p-value)"),
       colour = NULL,
       title = "Differential expression in 70 paired HCC samples",
       subtitle = paste0("DESeq2, paired design; FDR < ", padj_cutoff,
                         " and |log2FC| >= ", lfc_cutoff, "; top 10 genes each way labelled")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave("figures/03_volcano.png", volcano, width = 8, height = 7, dpi = 300)

# ---------- 4. MA plot ----------
ma <- ggplot(plot_df, aes(log10(baseMean), log2FC, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_colour_manual(values = dir_cols) +
  labs(x = expression(log[10]~"(mean normalised count)"),
       y = "Shrunken log2 fold change",
       colour = NULL,
       title = "MA plot: fold change vs expression level") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave("figures/04_ma_plot.png", ma, width = 8, height = 6, dpi = 300)

# ---------- Report: copy everything this prints and send it ----------
cat("\nCoefficients in the model (last one is the tumour effect):\n")
print(tail(resultsNames(dds), 3))
cat("\nGenes tested:", nrow(de), "| with an adjusted p-value:", sum(!is.na(de$padj)), "\n")
cat("FDR < 0.05 (any fold change):", sum(de$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR < 0.05 and |log2FC| >= 1): up", n_up, "| down", n_down, "\n")
cat("Genes with adjusted p-value too small to represent (shown capped):", sum(de$padj == 0, na.rm = TRUE), "\n")
show_cols <- c("symbol", "log2FC", "padj", "baseMean")
fmt <- function(x) { x$log2FC <- round(x$log2FC, 2); x$padj <- signif(x$padj, 2)
x$baseMean <- round(x$baseMean); x[, show_cols] }
cat("\nTop 15 UP in tumour (by significance):\n")
print(fmt(head(de[de$direction == "Up in tumour", ], 15)), row.names = FALSE)
cat("\nTop 15 DOWN in tumour (by significance):\n")
print(fmt(head(de[de$direction == "Down in tumour", ], 15)), row.names = FALSE)