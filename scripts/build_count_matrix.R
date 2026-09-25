#!/usr/bin/env Rscript
# Usage: Rscript build_count_matrix.R <data_dir>
# Reads all .genes.results files in <data_dir>/counts/ and builds two count matrices:
#   star_expected_count_ensembl.csv  (Ensembl IDs as row names)
#   star_expected_count_geneID.csv   (gene symbols as row names)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript build_count_matrix.R <data_dir>")

DATA_DIR <- args[1]
counts_dir <- file.path(DATA_DIR, "counts")

if (!dir.exists(counts_dir)) stop(paste("counts/ directory not found:", counts_dir))

suppressPackageStartupMessages(library(stringr))

files <- list.files(counts_dir)
if (length(files) == 0) stop("No .genes.results files found in counts/")

cat("Found", length(files), "samples\n")

# Initialize from first file
i <- files[1]
samplename <- strsplit(i, "[.]")[[1]][1]
data <- read.table(file.path(counts_dir, i), sep="\t", header=TRUE,
                   na.strings="", stringsAsFactors=FALSE)

gene_names_ensembl <- str_split_fixed(data$gene_id, "_", 2)[, 1]
star_counts_raw <- data.frame(row.names = gene_names_ensembl)

# Build matrix
for (i in files) {
  samplename <- strsplit(i, "[.]")[[1]][1]
  temp <- read.table(file.path(counts_dir, i), sep="\t", header=TRUE,
                     na.strings="", stringsAsFactors=FALSE)
  star_counts_raw[samplename] <- temp$expected_count
}

# Ensembl ID version
write.csv(star_counts_raw, file.path(DATA_DIR, "star_expected_count_ensembl.csv"))
save(star_counts_raw, file = file.path(DATA_DIR, "star_expected_count_ensembl.Rda"))
cat("Saved star_expected_count_ensembl.csv\n")

# Gene symbol version
gene_names_symbol <- str_split_fixed(data$gene_id, "_", 2)[, 2]
row.names(star_counts_raw) <- make.unique(gene_names_symbol, sep = ".")
write.csv(star_counts_raw, file.path(DATA_DIR, "star_expected_count_geneID.csv"))
save(star_counts_raw, file = file.path(DATA_DIR, "star_expected_count_geneID.Rda"))
cat("Saved star_expected_count_geneID.csv\n")
cat("Done. Count matrix dimensions:", nrow(star_counts_raw), "genes x", ncol(star_counts_raw), "samples\n")
