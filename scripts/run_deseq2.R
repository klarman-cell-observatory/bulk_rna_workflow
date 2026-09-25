#!/usr/bin/env Rscript
# Usage: Rscript run_deseq2.R <data_dir> <metadata_csv> <comparisons_json>
# comparisons_json: path to JSON file containing list of [group2, group1] pairs
#   group2 is numerator (condition of interest), group1 is reference/denominator
#   Output file naming: group2_vs_group1.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript run_deseq2.R <data_dir> <metadata_csv> <comparisons_json> [output_dir]")

DATA_DIR         <- args[1]
METADATA_PATH    <- args[2]
COMPARISONS_JSON <- args[3]
OUTPUT_DIR       <- if (length(args) >= 4) args[4] else file.path(DATA_DIR, "deseq2_results")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

# Load count matrix (gene symbol version)
count_file <- file.path(DATA_DIR, "star_expected_count_geneID.csv")
if (!file.exists(count_file)) stop(paste("Count matrix not found:", count_file))

cat("Loading count matrix...\n")
expression_matrix <- read.csv(count_file, stringsAsFactors=FALSE,
                               check.names=FALSE, row.names=1)

# Clean sample names (strip _S### suffix and replace . with -)
colnames(expression_matrix) <- gsub("_S[0-9]+.*$", "", colnames(expression_matrix))
colnames(expression_matrix) <- gsub("\\.", "-", colnames(expression_matrix))

# Load metadata
metadata <- read.csv(METADATA_PATH, stringsAsFactors=FALSE)
metadata$Sample_Name  <- factor(metadata$Sample_Name)
metadata$Sample_Group <- factor(metadata$Sample_Group)
rownames(metadata) <- metadata$Sample_Name

# Align columns to metadata
common_samples <- intersect(colnames(expression_matrix), rownames(metadata))
if (length(common_samples) == 0) stop("No samples match between count matrix and metadata")
cat("Samples matched:", length(common_samples), "\n")

expression_matrix <- expression_matrix[, common_samples]
metadata <- metadata[common_samples, ]

# Load comparisons
comparisons_raw <- readLines(COMPARISONS_JSON)
comparisons_json <- paste(comparisons_raw, collapse="")
# Parse simple JSON array of arrays: [[g2,g1],[g2,g1],...]
comparisons_json <- gsub("\\s+", "", comparisons_json)
inner <- gsub("^\\[\\[|\\]\\]$", "", comparisons_json)
pairs_str <- strsplit(inner, "\\],\\[")[[1]]
comparisons <- lapply(pairs_str, function(p) {
  vals <- gsub('"', '', strsplit(p, ",")[[1]])
  vals
})

cat("Running", length(comparisons), "comparisons...\n")
results_dir <- OUTPUT_DIR
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)

for (comp in comparisons) {
  group2 <- trimws(comp[1])  # numerator
  group1 <- trimws(comp[2])  # reference/denominator

  cat("\n---\nRunning:", group2, "vs", group1, "\n")

  selected <- metadata$Sample_Group %in% c(group1, group2)
  meta_sub <- metadata[selected, ]
  expr_sub <- expression_matrix[, rownames(meta_sub)]

  meta_sub$Sample_Group <- factor(meta_sub$Sample_Group, levels=c(group1, group2))
  cat("Samples - ", group1, ":", sum(meta_sub$Sample_Group == group1),
      " | ", group2, ":", sum(meta_sub$Sample_Group == group2), "\n")

  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData=round(expr_sub),
                                  colData=meta_sub,
                                  design=~Sample_Group)
    dds <- DESeq(dds)

    res       <- results(dds)
    res_shrink <- lfcShrink(dds, coef=2, type="apeglm")

    res        <- na.omit(res[order(res$log2FoldChange), ])
    res_shrink <- na.omit(res_shrink[order(res_shrink$log2FoldChange), ])

    comp_name <- paste0(group2, "_vs_", group1)
    write.csv(res,        file.path(results_dir, paste0(comp_name, ".csv")))
    write.csv(res_shrink, file.path(results_dir, paste0(comp_name, "_shrink.csv")))
    cat("Saved:", comp_name, "\n")
  }, error = function(e) {
    cat("ERROR for", group2, "vs", group1, ":", conditionMessage(e), "\n")
  })
}

cat("\nDESeq2 complete. Results written to:", results_dir, "\n")
