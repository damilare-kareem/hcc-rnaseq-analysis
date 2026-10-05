# 06_session_info.R - record R and package versions used for this analysis
suppressPackageStartupMessages({
  library(GEOquery); library(DESeq2); library(apeglm); library(ggplot2); library(ggrepel)
  library(pheatmap); library(clusterProfiler); library(enrichplot); library(org.Hs.eg.db); library(msigdbr)
})
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("Saved results/session_info.txt\n")