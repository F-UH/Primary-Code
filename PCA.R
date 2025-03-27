# ==============================
# Load Required Packages
# ==============================
library(dplyr)
library(stringr)
library(DESeq2)
library(limma)
library(ggplot2)

# ==============================
# Read Expression and Sample Data
# ==============================
# Load TPM expression data and sample annotations
datasets <- list(
  TCGA = list(tpm = "TCGA_tpm.csv", sample = "TCGA_sample.csv", type = "tissue"),
  GSE62182 = list(tpm = "GSE62182_tpm.csv", sample = "GSE62182_filter_sample.csv", type = "tissue"),
  GSE83527 = list(tpm = "GSE83527_tpm.csv", sample = "GSE83527_sample.csv", type = "tissue"),
  GSE175462 = list(tpm = "GSE175462_tpm.csv", sample = "GSE175462_sample.csv", type = "tissue"),
  GSE110907 = list(tpm = "GSE110907_tpm.csv", sample = "GSE110907_sample.csv", type = "tissue"),
  RUSH = list(tpm = "RUSH_tpm.csv", sample = "RUSH_sample.csv", type = "exosome")
)

expr_list <- list()
meta_list <- list()

for (name in names(datasets)) {
  expr <- read.csv(datasets[[name]]$tpm, header = TRUE, row.names = 1)
  sample <- read.csv(datasets[[name]]$sample, header = FALSE)
  sample$batch <- name
  sample$type <- datasets[[name]]$type
  
  # Filter out low-expression genes (50% or more samples have 0 TPM)
  expr <- expr[rowSums(expr == 0) / ncol(expr) < 0.5, ]
  
  expr_list[[name]] <- expr
  meta_list[[name]] <- sample
}

# ==============================
# Identify Common Genes
# ==============================
# Get the intersection of genes across all datasets
common_genes <- Reduce(intersect, lapply(expr_list, rownames))

# Filter to common genes only
expr_list <- lapply(expr_list, function(x) x[common_genes, ])

# Combine all expression and sample data
all_expr <- do.call(cbind, expr_list)
all_meta <- do.call(rbind, meta_list)
all_expr_matrix <- t(all_expr)

# ==============================
# Batch Effect Removal with limma
# ==============================
# Design matrix using sample diagnosis (assumed to be in V2)
design <- model.matrix(~ V2, data = all_meta)

# Remove batch effect
expr_no_batch <- removeBatchEffect(all_expr, batch = all_meta$batch, design = design)

# ==============================
# PCA Analysis
# ==============================
# Prepare PCA input
expr_scaled <- scale(t(expr_no_batch))  # transpose to samples x features
pc <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)

# Check if any NAs introduced
stopifnot(!any(is.na(pc$x)))

# Create PCA result dataframe
pc_results <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                         group = all_meta$V2,
                         batch = factor(all_meta$batch, 
                                        levels = c("TCGA", "GSE83527", "GSE62182", "GSE175462", "GSE110907", "RUSH")))

# ==============================
# Plot PCA
# ==============================
svg("PCA_all_limma.svg", width = 7, height = 5)
pca_plot <- ggplot(pc_results, aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 1.5) +
  stat_ellipse(aes(group = group), geom = "polygon", alpha = 0.2, fill = NA) +
  xlab("PC 1") + ylab("PC 2") +
  labs(title = "PCA Analysis of tRNA Datasets",
       shape = "Batch", color = "Diagnosis") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black"),
    legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid"),
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(19, 15, 17, 18, 20, 13)) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(-20, 40))
print(pca_plot)
dev.off()
