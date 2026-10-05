# 05_enrichment.R
# Which biological processes and pathways do the differentially expressed genes point to?
#  1. GO biological processes (over-representation, up and down genes separately)
#  2. KEGG pathways (over-representation)
#  3. MSigDB Hallmark gene sets (GSEA on all genes, ranked)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(ggplot2)
})

de <- readRDS("data/processed/de_results.rds")

up_genes   <- unique(de$ensembl[de$direction == "Up in tumour"])
down_genes <- unique(de$ensembl[de$direction == "Down in tumour"])
universe   <- unique(de$ensembl)   # every gene that was tested = the fair background
gene_groups <- list(`Up in tumour` = up_genes, `Down in tumour` = down_genes)

# ---------- 1. GO biological processes ----------
cat("GO enrichment, started", format(Sys.time(), "%H:%M"), "\n")
go <- compareCluster(geneClusters = gene_groups, fun = "enrichGO",
                     universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
go <- simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)   # merge near-duplicate terms
go_df <- as.data.frame(go)
write.csv(go_df, "results/go_bp_enrichment.csv", row.names = FALSE)

p_go <- dotplot(go, showCategory = 10, label_format = 45) +
  labs(title = "GO biological processes", subtitle = "Top 10 terms per group (FDR < 0.05)")
ggsave("figures/08_go_enrichment.png", p_go, width = 9, height = 9, dpi = 300)

# ---------- 2. KEGG pathways (needs internet: data comes from the KEGG website) ----------
cat("KEGG enrichment, started", format(Sys.time(), "%H:%M"), "\n")
to_entrez <- function(ids) {
  unique(na.omit(mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID",
                        keytype = "ENSEMBL", multiVals = "first")))
}
kegg <- tryCatch(
  compareCluster(geneClusters = lapply(gene_groups, to_entrez), fun = "enrichKEGG",
                 organism = "hsa", universe = to_entrez(universe), pvalueCutoff = 0.05),
  error = function(e) { message("KEGG step skipped: ", conditionMessage(e)); NULL })

kegg_df <- NULL
if (!is.null(kegg)) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg)
  write.csv(kegg_df, "results/kegg_enrichment.csv", row.names = FALSE)
  p_kegg <- dotplot(kegg, showCategory = 10, label_format = 45) +
    labs(title = "KEGG pathways", subtitle = "Top 10 pathways per group (FDR < 0.05)")
  ggsave("figures/09_kegg_enrichment.png", p_kegg, width = 9, height = 8, dpi = 300)
}

# ---------- 3. Hallmark GSEA ----------
cat("Hallmark GSEA, started", format(Sys.time(), "%H:%M"), "\n")
hallmark <- tryCatch(msigdbr(collection = "H"),                 # current msigdbr
                     error = function(e) msigdbr(category = "H")) # older versions
t2g <- as.data.frame(unique(hallmark[, c("gs_name", "gene_symbol")]))

# Rank every tested gene by its Wald statistic (positive = higher in tumour)
ranked <- de[!is.na(de$stat), c("symbol", "stat")]
ranked <- ranked[order(-abs(ranked$stat)), ]
ranked <- ranked[!duplicated(ranked$symbol), ]       # one value per gene symbol
gene_list <- sort(setNames(ranked$stat, ranked$symbol), decreasing = TRUE)

set.seed(42)
gsea <- GSEA(gene_list, TERM2GENE = t2g, pvalueCutoff = 0.05,
             eps = 0, seed = TRUE, verbose = FALSE)
hm <- as.data.frame(gsea)
write.csv(hm, "results/hallmark_gsea.csv", row.names = FALSE)

hm$pathway   <- gsub("_", " ", sub("^HALLMARK_", "", hm$ID))
hm$direction <- ifelse(hm$NES > 0, "Higher in tumour", "Lower in tumour")
p_hm <- ggplot(hm, aes(NES, reorder(pathway, NES), fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c("Higher in tumour" = "#D55E00", "Lower in tumour" = "#0072B2")) +
  labs(x = "Normalised enrichment score (NES)", y = NULL, fill = NULL,
       title = "Hallmark gene sets in HCC vs matched non-tumour liver",
       subtitle = paste0("GSEA on ", length(gene_list),
                         " genes ranked by DESeq2 Wald statistic; FDR < 0.05")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/10_hallmark_gsea.png", p_hm, width = 9, height = 2 + 0.22 * nrow(hm), dpi = 300)

# Running-score plot for the strongest set in each direction
top_pos <- hm$ID[which.max(hm$NES)]
top_neg <- hm$ID[which.min(hm$NES)]
p_run <- gseaplot2(gsea, geneSetID = c(top_pos, top_neg), color = c("#D55E00", "#0072B2"),
                   title = "Strongest Hallmark gene sets in each direction")
png("figures/11_gsea_running_score.png", width = 9, height = 6, units = "in", res = 300)
print(p_run)
dev.off()

# ---------- Report: copy everything this prints and send it ----------
options(width = 150)
top_terms <- function(df, group, n = 10) {
  x <- df[df$Cluster == group, ]
  x <- x[order(x$p.adjust), c("Description", "Count", "p.adjust")]
  x$p.adjust <- signif(x$p.adjust, 2)
  head(x, n)
}
cat("\nGenes used: up", length(up_genes), "| down", length(down_genes), "| background", length(universe), "\n")
cat("\nGO terms (after merging near-duplicates): up", sum(go_df$Cluster == "Up in tumour"),
    "| down", sum(go_df$Cluster == "Down in tumour"), "\n")
cat("\nTop GO terms, UP in tumour:\n");   print(top_terms(go_df, "Up in tumour"), row.names = FALSE)
cat("\nTop GO terms, DOWN in tumour:\n"); print(top_terms(go_df, "Down in tumour"), row.names = FALSE)
if (!is.null(kegg_df)) {
  cat("\nTop KEGG pathways, UP in tumour:\n");   print(top_terms(kegg_df, "Up in tumour"), row.names = FALSE)
  cat("\nTop KEGG pathways, DOWN in tumour:\n"); print(top_terms(kegg_df, "Down in tumour"), row.names = FALSE)
}
cat("\nHallmark sets significant:", nrow(hm), "\n")
hm_show <- hm[order(-hm$NES), c("pathway", "setSize", "NES", "p.adjust")]
hm_show$NES <- round(hm_show$NES, 2); hm_show$p.adjust <- signif(hm_show$p.adjust, 2)
print(hm_show, row.names = FALSE)