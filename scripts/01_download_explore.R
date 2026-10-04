# 01_download_explore.R
# Download GSE144269 (70 paired HCC tumour / non-tumour liver RNA-seq) and inspect it

library(GEOquery)
options(timeout = 600)   # allow slow downloads

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 1. Supplementary files: the processed gene count tables
supp <- getGEOSuppFiles("GSE144269", baseDir = "data/raw")
files <- rownames(supp)

# 2. Sample information (tissue type, patient ID, etc.)
gse <- getGEO("GSE144269", GSEMatrix = TRUE, destdir = "data/raw")
pheno <- pData(gse[[1]])
saveRDS(pheno, "data/raw/pheno.rds")

# 3. Read the RSEM gene count table
counts_file <- files[grepl("count", basename(files), ignore.case = TRUE)][1]
raw_counts <- read.delim(counts_file, check.names = FALSE)
saveRDS(raw_counts, "data/raw/raw_counts.rds")

# ---------- Report: copy everything below this line's output and send it ----------
cat(R.version.string, "\n\n")
cat("Supplementary files:\n"); print(basename(files))
cat("\nCount table read from:", basename(counts_file), "\n")
cat("Dimensions (rows x columns):", dim(raw_counts), "\n\n")
print(raw_counts[1:5, 1:4])
cat("\nSample metadata:", nrow(pheno), "samples\n")
print(head(pheno[, c("title", "geo_accession", "source_name_ch1")], 6))
ch1 <- grep(":ch1$", colnames(pheno), value = TRUE)
for (field in ch1) cat("\n", field, ":", paste(head(unique(pheno[[field]]), 5), collapse = " | "))