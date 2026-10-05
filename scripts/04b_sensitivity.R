# 04b_sensitivity.R
# Sensitivity check: are the results driven by a few low-quality samples?
# Pre-specified rule: remove every patient who has a library below 5 million reads
# (both of their samples, to keep the design paired), refit DESeq2, and compare.

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(ggplot2)
})

n_workers <- min(4, max(1, parallel::detectCores() - 1))
bp <- if (n_workers > 1) SnowParam(workers = n_workers) else SerialParam()

dds_all    <- readRDS("data/processed/dds.rds")          # filtered, not yet fitted
de_full    <- readRDS("data/processed/de_results.rds")
raw_counts <- readRDS("data/raw/raw_counts.rds")

# ---------- 1. Apply the rule ----------
min_reads <- 5e6
lib_total <- colSums(raw_counts[, -1])                   # total reads per sample, before filtering
low_samples   <- names(lib_total)[lib_total < min_reads]
drop_patients <- unique(as.character(dds_all$patient[colnames(dds_all) %in% low_samples]))

dds_red <- dds_all[, !(dds_all$patient %in% drop_patients)]
dds_red$patient <- droplevels(dds_red$patient)

# ---------- 2. Refit (slow; saved so a re-run skips it) ----------
fitted_file <- "data/processed/dds_sensitivity_fitted.rds"
if (file.exists(fitted_file)) {
  dds_red <- readRDS(fitted_file)
  cat("Loaded the saved sensitivity model.\n")
} else {
  cat("Fitting DESeq2 without", length(drop_patients), "patient(s), started",
      format(Sys.time(), "%H:%M"), "\n")
  dds_red <- DESeq(dds_red, parallel = n_workers > 1, BPPARAM = bp)
  saveRDS(dds_red, fitted_file)
  cat("Model fitted at", format(Sys.time(), "%H:%M"), "\n")
}
res_red <- results(dds_red, name = "tissue_tumor_vs_non_tumor", alpha = 0.05)

# ---------- 3. Compare with the full analysis ----------
# Both sides use unshrunk fold changes and the same thresholds, so the comparison is like for like
cmp <- merge(de_full[, c("gene_key", "symbol", "log2FC_unshrunk", "padj")],
             data.frame(gene_key = rownames(res_red),
                        log2FC_red = res_red$log2FoldChange,
                        padj_red = res_red$padj),
             by = "gene_key")
is_sig <- function(lfc, p) !is.na(p) & p < 0.05 & abs(lfc) >= 1
cmp$sig_full <- is_sig(cmp$log2FC_unshrunk, cmp$padj)
cmp$sig_red  <- is_sig(cmp$log2FC_red, cmp$padj_red)

r_lfc    <- cor(cmp$log2FC_unshrunk, cmp$log2FC_red, use = "complete.obs")
n_full   <- sum(cmp$sig_full)
n_red    <- sum(cmp$sig_red)
n_both   <- sum(cmp$sig_full & cmp$sig_red)
jaccard  <- n_both / sum(cmp$sig_full | cmp$sig_red)
same_dir <- mean(sign(cmp$log2FC_unshrunk[cmp$sig_full]) == sign(cmp$log2FC_red[cmp$sig_full]))

top50 <- head(de_full$gene_key[de_full$direction != "Not significant"], 50)
top50_kept <- sum(cmp$sig_red[match(top50, cmp$gene_key)])

marker_symbols <- c("GPC3", "AKR1B10", "AFP", "MKI67",
                    "VIPR1", "CYP1A2", "FCN3", "ECM1", "LIFR", "STAB2", "CLEC4G")
marker_cmp <- cmp[match(marker_symbols, cmp$symbol), c("symbol", "log2FC_unshrunk", "padj", "log2FC_red", "padj_red")]
marker_cmp <- marker_cmp[!is.na(marker_cmp$symbol), ]
marker_cmp[, c("log2FC_unshrunk", "log2FC_red")] <- round(marker_cmp[, c("log2FC_unshrunk", "log2FC_red")], 2)
marker_cmp[, c("padj", "padj_red")] <- signif(marker_cmp[, c("padj", "padj_red")], 2)

write.csv(cmp, "results/sensitivity_comparison.csv", row.names = FALSE)

# ---------- 4. Figure ----------
sens_plot <- ggplot(cmp, aes(log2FC_unshrunk, log2FC_red)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 0.6, alpha = 0.3, colour = "grey30") +
  coord_equal() +
  labs(x = paste0("log2 fold change, all ", nlevels(dds_all$patient), " patients"),
       y = paste0("log2 fold change, ", nlevels(dds_red$patient), " patients"),
       title = "Sensitivity analysis: excluding low-depth samples",
       subtitle = paste0("Removed ", length(drop_patients), " patient(s) with a library < 5M reads\n",
                         "Pearson r = ", formatC(r_lfc, digits = 3, format = "f"),
                         " across ", nrow(cmp), " genes")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave("figures/07_sensitivity.png", sens_plot, width = 7, height = 7, dpi = 300)

# ---------- Report: copy everything this prints and send it ----------
options(width = 150)
cat("\nSamples below 5M reads:\n"); print(round(sort(lib_total[low_samples]) / 1e6, 2))
cat("Patients removed (both samples each):", paste(drop_patients, collapse = ", "), "\n")
cat("Patients remaining:", nlevels(dds_red$patient), "\n\n")
cat("Correlation of fold changes (all genes): r =", round(r_lfc, 4), "\n")
cat("Significant genes (FDR < 0.05, |unshrunk log2FC| >= 1): full", n_full, "| reduced", n_red,
    "| in both", n_both, "\n")
cat("Full-analysis significant genes still significant:", round(100 * n_both / n_full, 1), "%\n")
cat("Jaccard overlap of the two significant sets:", round(jaccard, 3), "\n")
cat("Same direction of change (full-analysis significant genes):", round(100 * same_dir, 2), "%\n")
cat("Top 50 genes still significant:", top50_kept, "of 50\n")
cat("\nMarker genes, full vs reduced:\n"); print(marker_cmp, row.names = FALSE)