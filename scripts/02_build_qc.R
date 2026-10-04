# 02_build_qc.R
# Build the sample table and count matrix, pair tumours with their normal liver,
# filter low-count genes, and run quality control (PCA, sample clustering).

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dir.create("data/processed", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

raw_counts <- readRDS("data/raw/raw_counts.rds")
pheno      <- readRDS("data/raw/pheno.rds")

# ---------- 1. Sample table ----------
# Titles look like "pat 10 tumor [HDV_10A_S33]": patient number + tissue + sample ID
sample_table <- data.frame(
  sample_id = sub(".*\\[(.*)\\].*", "\\1", pheno$title),
  gsm       = pheno$geo_accession,
  patient   = paste0("P", sub("^pat (\\d+) .*", "\\1", pheno$title)),
  tissue    = ifelse(pheno[["tumor/non-tumor:ch1"]] == "tumor", "tumor", "non_tumor")
)
sample_table$tissue  <- factor(sample_table$tissue, levels = c("non_tumor", "tumor"))
sample_table$patient <- factor(sample_table$patient)
sample_table <- sample_table[order(sample_table$patient, sample_table$tissue), ]
rownames(sample_table) <- sample_table$sample_id

# ---------- 2. Integrity checks ----------
check <- function(ok, msg) {
  cat(if (ok) "PASS" else "FAIL", "-", msg, "\n")
  ok
}
count_ids  <- colnames(raw_counts)[-1]
id_patient <- paste0("P", sub("^HDV_(\\d+)[AB]_.*", "\\1", sample_table$sample_id))
id_letter  <- sub("^HDV_\\d+([AB])_.*", "\\1", sample_table$sample_id)

checks <- c(
  check(setequal(count_ids, sample_table$sample_id),
        "every count-table column matches exactly one sample in the metadata"),
  check(!anyDuplicated(sample_table$sample_id), "sample IDs are unique"),
  check(all(id_patient == as.character(sample_table$patient)),
        "patient number in sample ID agrees with the title"),
  check(all((id_letter == "A") == (sample_table$tissue == "tumor")),
        "A = tumour and B = non-tumour in every sample ID"),
  check(all(table(sample_table$patient) == 2) &&
          all(table(sample_table$patient, sample_table$tissue) == 1),
        "every patient has exactly one tumour and one non-tumour sample")
)
if (!all(checks)) stop("Integrity check failed - send the output above.")

# ---------- 3. Count matrix + gene annotation ----------
# Gene IDs look like "ENSG00000000003.14|TSPAN6"
gene_key <- sub("\\|.*$", "", raw_counts[[1]])
genes <- data.frame(
  gene_key = make.unique(gene_key),
  ensembl  = sub("\\..*$", "", gene_key),
  symbol   = sub("^[^|]*\\|", "", raw_counts[[1]])
)

counts <- as.matrix(raw_counts[, -1])
rownames(counts) <- genes$gene_key
counts <- counts[, sample_table$sample_id]   # same order as the sample table
non_integer <- mean(counts != round(counts))
counts <- round(counts)                      # RSEM gives expected counts; DESeq2 needs integers
mode(counts) <- "integer"

# ---------- 4. DESeq2 object with a paired design ----------
# "~ patient + tissue": compare tumour with normal liver *within* each patient
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_table,
                              design = ~ patient + tissue)
rowData(dds)$ensembl <- genes$ensembl
rowData(dds)$symbol  <- genes$symbol

library_size_m <- colSums(counts(dds)) / 1e6
genes_before   <- nrow(dds)

# Keep genes with at least 10 reads in at least 70 samples (the size of the smallest group)
keep <- rowSums(counts(dds) >= 10) >= 70
dds  <- dds[keep, ]

# ---------- 5. Quality control ----------
vsd <- vst(dds, blind = TRUE)   # variance-stabilised values, for plots only

tissue_cols   <- c(non_tumor = "#0072B2", tumor = "#D55E00")
tissue_labels <- c(non_tumor = "Non-tumour liver", tumor = "HCC tumour")

# PCA, with a grey line joining each patient's two samples
pca <- plotPCA(vsd, intgroup = c("tissue", "patient"), ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
pca_plot <- ggplot(pca, aes(PC1, PC2)) +
  geom_line(aes(group = patient), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = tissue), size = 2.5, alpha = 0.9) +
  scale_colour_manual(values = tissue_cols, labels = tissue_labels) +
  labs(x = paste0("PC1 (", pv[1], "% variance)"),
       y = paste0("PC2 (", pv[2], "% variance)"),
       colour = NULL,
       title = "PCA of 140 liver samples",
       subtitle = "Top 500 most variable genes; grey lines join each patient's tumour and non-tumour sample") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/01_pca.png", pca_plot, width = 8, height = 6, dpi = 300)

# Sample-to-sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
annotation  <- data.frame(Tissue = tissue_labels[as.character(vsd$tissue)],
                          row.names = colnames(vsd))
pheatmap(as.matrix(sample_dist),
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = annotation, annotation_row = annotation,
         annotation_colors = list(Tissue = setNames(tissue_cols, tissue_labels)),
         show_rownames = FALSE, show_colnames = FALSE,
         color = hcl.colors(100, "Blues 3"),
         main = "Sample-to-sample distances (140 samples)",
         filename = "figures/02_sample_distances.png", width = 8, height = 7)

# How well does unsupervised clustering recover tumour vs non-tumour?
clusters <- cutree(hclust(sample_dist), k = 2)

# ---------- 6. Save ----------
saveRDS(dds, "data/processed/dds.rds")
saveRDS(vsd, "data/processed/vsd.rds")
write.csv(sample_table, "results/sample_table.csv", row.names = FALSE)

# ---------- Report: copy everything this prints and send it ----------
cat("\nSamples per tissue:\n"); print(table(sample_table$tissue))
cat("Patients with a complete pair:", nlevels(sample_table$patient), "\n")
cat("Fraction of non-integer counts before rounding:", round(non_integer, 4), "\n")
cat("Library size (million reads): min", round(min(library_size_m), 1),
    "| median", round(median(library_size_m), 1), "| max", round(max(library_size_m), 1), "\n")
cat("Smallest 3 libraries:\n"); print(round(sort(library_size_m)[1:3], 2))
cat("Genes before filtering:", genes_before, "| after filtering:", nrow(dds), "\n")
cat("PCA variance explained: PC1", pv[1], "% | PC2", pv[2], "%\n")
cat("\nUnsupervised clustering (2 clusters) vs tissue:\n")
print(table(cluster = clusters, tissue = vsd$tissue))
