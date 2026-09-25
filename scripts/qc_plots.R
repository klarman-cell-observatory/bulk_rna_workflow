#!/usr/bin/env Rscript
# Usage: Rscript qc_plots.R <data_dir> <metadata_csv>
# Produces:
#   Figures/pca_plot.png
#   Figures/correlation_heatmap.png
#   Figures/qc_flags.csv   (machine-readable flags for the notebook to display)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript qc_plots.R <data_dir> <metadata_csv>")

DATA_DIR      <- args[1]
METADATA_PATH <- args[2]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

fig_dir <- file.path(DATA_DIR, "Figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading count matrix...\n")
count_file <- file.path(DATA_DIR, "star_expected_count_geneID.csv")
expr <- read.csv(count_file, row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)
colnames(expr) <- gsub("_S[0-9]+.*$", "", colnames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))

meta <- read.csv(METADATA_PATH, stringsAsFactors = FALSE)
rownames(meta) <- meta$Sample_Name
meta$Sample_Group <- factor(meta$Sample_Group)

common <- intersect(colnames(expr), rownames(meta))
if (length(common) == 0) stop("No samples match between count matrix and metadata")
expr <- expr[, common]
meta <- meta[common, ]
cat("Samples:", length(common), "\n")

# ── VST normalization ─────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = round(expr),
                              colData   = meta,
                              design    = ~ Sample_Group)
dds <- estimateSizeFactors(dds)

# Flag samples with extreme size factors
sf     <- sizeFactors(dds)
sf_med <- median(sf)
sf_flags <- names(sf)[sf > sf_med * 2 | sf < sf_med / 2]

vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

# ── PCA ───────────────────────────────────────────────────────────────────────
cat("Computing PCA...\n")
pca     <- prcomp(t(mat), scale. = FALSE)
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- data.frame(
  PC1          = pca$x[, 1],
  PC2          = pca$x[, 2],
  Sample       = rownames(pca$x),
  Sample_Group = meta[rownames(pca$x), "Sample_Group"]
)

# Detect PCA outliers: samples > 2.5 SD from their group centroid on PC1+PC2
pca_flags <- character(0)
for (grp in levels(pca_df$Sample_Group)) {
  grp_data <- pca_df[pca_df$Sample_Group == grp, ]
  if (nrow(grp_data) < 3) next
  for (pc in c("PC1", "PC2")) {
    mu  <- mean(grp_data[[pc]])
    ssd <- sd(grp_data[[pc]])
    if (is.na(ssd) || ssd == 0) next
    outliers <- grp_data$Sample[abs(grp_data[[pc]] - mu) > 2.5 * ssd]
    pca_flags <- union(pca_flags, outliers)
  }
}

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample_Group, label = Sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
  labs(
    x     = paste0("PC1 (", pct_var[1], "% variance)"),
    y     = paste0("PC2 (", pct_var[2], "% variance)"),
    title = "PCA — VST-normalized counts",
    subtitle = "Replicates should cluster together; conditions should separate"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# Add red circles around flagged samples
if (length(pca_flags) > 0) {
  flag_df <- pca_df[pca_df$Sample %in% pca_flags, ]
  p <- p + geom_point(data = flag_df, size = 6, shape = 1, color = "red", stroke = 1.5)
}

ggsave(file.path(fig_dir, "pca_plot.png"), p, width = 9, height = 6, dpi = 200)
cat("Saved pca_plot.png\n")

# ── Sample correlation heatmap ────────────────────────────────────────────────
cat("Computing sample correlations...\n")
cor_mat  <- cor(mat, method = "pearson")
cor_flag_thresh <- 0.90

# Flag samples whose minimum intra-group correlation is below threshold
cor_flags <- character(0)
for (grp in levels(meta$Sample_Group)) {
  grp_samps <- rownames(meta)[meta$Sample_Group == grp]
  grp_samps <- intersect(grp_samps, colnames(cor_mat))
  if (length(grp_samps) < 2) next
  sub <- cor_mat[grp_samps, grp_samps]
  diag(sub) <- NA
  min_cors  <- apply(sub, 1, min, na.rm = TRUE)
  bad       <- names(min_cors)[min_cors < cor_flag_thresh]
  cor_flags <- union(cor_flags, bad)
}

# Annotation sidebar
ann_col <- data.frame(Group = meta$Sample_Group)
rownames(ann_col) <- rownames(meta)

n_groups <- nlevels(meta$Sample_Group)
pal      <- colorRampPalette(brewer.pal(min(n_groups, 8), "Set2"))(n_groups)
ann_colors <- list(Group = setNames(pal, levels(meta$Sample_Group)))

heat_colors <- colorRampPalette(c("#2166ac", "#f7f7f7", "#d6604d"))(100)

png(file.path(fig_dir, "correlation_heatmap.png"),
    width = 1400, height = 1200, res = 150)
pheatmap(cor_mat,
         color            = heat_colors,
         #breaks           = seq(0.7, 1.0, length.out = 101),
         annotation_col   = ann_col,
         annotation_colors = ann_colors,
         show_rownames    = TRUE,
         show_colnames    = TRUE,
         fontsize         = 7,
         main             = "Sample-to-sample Pearson correlation (VST counts)")
dev.off()
cat("Saved correlation_heatmap.png\n")

# ── Write flags CSV ───────────────────────────────────────────────────────────
sf_detail  <- if (length(sf_flags)  > 0) paste0("Size factor ", round(sf[sf_flags], 2), " vs median ", round(sf_med, 2)) else character(0)
pca_detail <- rep("Outlier (>2.5 SD from group centroid on PC1 or PC2)", length(pca_flags))
cor_detail <- rep(paste0("Min intra-group Pearson r < ", cor_flag_thresh), length(cor_flags))

all_flags <- data.frame(
  sample  = c(sf_flags,  pca_flags,  cor_flags),
  check   = c(rep("size_factor",               length(sf_flags)),
              rep("pca_outlier",               length(pca_flags)),
              rep("low_intragroup_correlation", length(cor_flags))),
  detail  = c(sf_detail, pca_detail, cor_detail),
  stringsAsFactors = FALSE
)

write.csv(all_flags, file.path(fig_dir, "qc_flags.csv"), row.names = FALSE)
cat("QC flags:", nrow(all_flags), "\n")
cat("Done.\n")
