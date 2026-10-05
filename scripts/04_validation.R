# 04_validation.R
# Does the pipeline recover known HCC biology?
#  1. Heatmap of the top 50 differentially expressed genes
#  2. Published HCC marker genes: does each move in the expected direction?
#  3. Identify the most unusual samples (for the sensitivity check in Bit 4b)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dds <- readRDS("data/processed/dds_fitted.rds")
vsd <- readRDS("data/processed/vsd.rds")
de  <- readRDS("data/processed/de_results.rds")

tissue_cols   <- c(non_tumor = "#0072B2", tumor = "#D55E00")
tissue_labels <- c(non_tumor = "Non-tumour liver", tumor = "HCC tumour")

# ---------- 1. Heatmap of the top 50 genes ----------
top_up   <- head(de$gene_key[de$direction == "Up in tumour"], 25)
top_down <- head(de$gene_key[de$direction == "Down in tumour"], 25)
top_genes <- c(top_up, top_down)

z <- assay(vsd)[top_genes, ]
z <- t(scale(t(z)))            # z-score each gene across samples
z[z > 3] <- 3; z[z < -3] <- -3 # cap extreme values so the colours stay readable
rownames(z) <- de$symbol[match(top_genes, de$gene_key)]

# Order columns: all non-tumour samples, then all tumours; cluster within each group
col_order <- unlist(lapply(c("non_tumor", "tumor"), function(t) {
  s <- colnames(vsd)[vsd$tissue == t]
  s[hclust(dist(t(z[, s])))$order]
}))
annotation <- data.frame(Tissue = tissue_labels[as.character(vsd$tissue)],
                         row.names = colnames(vsd))
gene_annotation <- data.frame(Direction = rep(c("Up in tumour", "Down in tumour"), each = 25),
                              row.names = rownames(z))

pheatmap(z[, col_order],
         cluster_cols = FALSE, cluster_rows = TRUE,
         gaps_col = sum(vsd$tissue == "non_tumor"),
         annotation_col = annotation, annotation_row = gene_annotation,
         annotation_colors = list(
           Tissue    = setNames(tissue_cols, tissue_labels),
           Direction = c("Up in tumour" = "#D55E00", "Down in tumour" = "#0072B2")),
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102),
         show_colnames = FALSE, fontsize_row = 7, border_color = NA,
         main = "Top 50 differentially expressed genes (z-scored expression)",
         filename = "figures/05_top50_heatmap.png", width = 10, height = 9)

# ---------- 2. Published HCC markers ----------
# Chosen from the literature, with the direction reported there:
#   GPC3: absent in normal adult liver, high in HCC (Haruyama & Kataoka 2016)
#   AKR1B10: higher in tumour than peri-tumour liver (Sci Rep 2016)
#   AFP: classic HCC serum marker, raised in a subset of tumours
#   MKI67: Ki-67, the standard proliferation marker
#   VIPR1, CYP1A2, FCN3, ECM1, LIFR: consistently low in HCC (Ogholbake & Cheng 2025)
#   STAB2, CLEC4G: sinusoidal endothelial markers lost in chronic liver disease (Verhulst et al. 2021)
markers <- data.frame(
  symbol   = c("GPC3", "AKR1B10", "AFP", "MKI67",
               "VIPR1", "CYP1A2", "FCN3", "ECM1", "LIFR", "STAB2", "CLEC4G"),
  expected = c(rep("Up", 4), rep("Down", 7))
)

norm_counts <- counts(dds, normalized = TRUE)
patients <- levels(dds$patient)
tumour_col <- sapply(patients, function(p) colnames(dds)[dds$patient == p & dds$tissue == "tumor"])
normal_col <- sapply(patients, function(p) colnames(dds)[dds$patient == p & dds$tissue == "non_tumor"])

paired_rows <- list()
for (i in seq_len(nrow(markers))) {
  g   <- markers$symbol[i]
  hit <- match(g, de$symbol)
  if (is.na(hit)) {
    markers$log2FC[i] <- NA; markers$padj[i] <- NA; markers$pct_patients[i] <- NA
    markers$verdict[i] <- "Not tested (filtered: too few reads)"
    next
  }
  key <- de$gene_key[hit]
  markers$log2FC[i] <- round(de$log2FC[hit], 2)
  markers$padj[i]   <- signif(de$padj[hit], 2)
  
  # Per-patient log2(tumour / matched normal)
  ratio <- log2((norm_counts[key, tumour_col] + 1) / (norm_counts[key, normal_col] + 1))
  in_expected <- if (markers$expected[i] == "Up") ratio > 0 else ratio < 0
  markers$pct_patients[i] <- round(100 * mean(in_expected))
  paired_rows[[g]] <- data.frame(symbol = g, expected = markers$expected[i], ratio = ratio)
  
  significant <- !is.na(de$padj[hit]) && de$padj[hit] < 0.05
  right_way   <- sign(de$log2FC[hit]) == ifelse(markers$expected[i] == "Up", 1, -1)
  markers$verdict[i] <- if (significant && right_way) "Agrees" else
    if (significant) "Opposite direction" else "Not significant"
}
write.csv(markers, "results/marker_validation.csv", row.names = FALSE)

# Figure: each dot is one patient's tumour vs own normal liver
paired <- do.call(rbind, paired_rows)
paired$symbol   <- factor(paired$symbol, levels = markers$symbol)
paired$expected <- factor(paste("Expected", tolower(paired$expected), "in HCC"),
                          levels = c("Expected up in HCC", "Expected down in HCC"))
marker_plot <- ggplot(paired, aes(symbol, ratio)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.2, height = 0, size = 1.2, alpha = 0.5, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "black", width = 0.5, linewidth = 0.4) +
  facet_grid(. ~ expected, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "log2 (tumour / matched non-tumour liver)",
       title = "Published HCC marker genes in 70 patients",
       subtitle = "Each dot is one patient: tumour expression relative to their own non-tumour liver") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(face = "italic"))
ggsave("figures/06_marker_validation.png", marker_plot, width = 10, height = 5.5, dpi = 300)

# ---------- 3. Most unusual samples ----------
d <- as.matrix(dist(t(assay(vsd))))
med_dist <- apply(d, 1, function(x) median(x[x > 0]))
depth <- colSums(counts(dds)) / 1e6
top_odd <- names(sort(med_dist, decreasing = TRUE))[1:5]
unusual <- data.frame(sample = top_odd,
                      tissue = vsd$tissue[match(top_odd, colnames(vsd))],
                      median_distance = round(med_dist[top_odd]),
                      million_reads = round(depth[top_odd], 1),
                      row.names = NULL)

# ---------- Report: copy everything this prints and send it ----------
options(width = 150)
cat("\nPublished marker check:\n")
print(markers, row.names = FALSE)
cat("\nMarkers agreeing with the literature:", sum(markers$verdict == "Agrees"),
    "of", sum(!grepl("Not tested", markers$verdict)), "tested\n")
cat("\nTypical sample's median distance to the others:", round(median(med_dist)), "\n")
cat("Five most unusual samples:\n"); print(unusual)