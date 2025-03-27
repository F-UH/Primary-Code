# ===============================
# Load Required Packages
# ===============================
library(pheatmap)
library(RColorBrewer)
library(corrplot)
library(limma)

# ===============================
# Define Plotting Function
# ===============================
plot_correlation_heatmap <- function(cor_matrix, filename, title) {
  colors <- colorRampPalette(brewer.pal(15, "RdBu"))(255)
  svg(filename, width = 10, height = 8)
  pheatmap(cor_matrix,
           col = colors,
           fontsize_row = 10,
           fontsize_col = 10,
           main = title,
           cluster_rows = FALSE, cluster_cols = FALSE,
           cellwidth = 20, cellheight = 20,
           legend = TRUE,
           breaks = seq(-1, 1, length.out = length(colors)),
           show_rownames = TRUE, show_colnames = TRUE)
  dev.off()
}

plot_corrplot <- function(cor_matrix, filename, title = "") {
  colors <- colorRampPalette(brewer.pal(10, "RdYlBu"))(255)
  svg(filename, width = 10, height = 8)
  corrplot(cor_matrix,
           method = "shade",
           type = "upper",
           addCoef.col = "white",
           number.cex = 0.7,
           col = colors,
           tl.cex = 0.8,
           tl.col = "black",
           tl.srt = 45)
  dev.off()
}

# ===============================
# Dataset-wise Spearman Correlations
# ===============================
datasets <- list(
  RUSH = "RUSH_filter.csv",
  TCGA = "TCGA_filter.csv",
  GSE110907 = "GSE110907_filter.csv",
  GSE62182 = "GSE62182_filter.csv",
  GSE83527 = "GSE83527_filter.csv",
  GSE175462 = "GSE175462_filter.csv"
)

for (name in names(datasets)) {
  data <- read.csv(datasets[[name]], header = TRUE, row.names = 1)
  data_t <- t(data)
  spearman_matrix <- cor(data_t, method = "spearman")
  
  # Heatmap
  plot_correlation_heatmap(spearman_matrix, paste0(name, ".svg"), 
                           paste("Gene Expression Correlation of", name))
  
  # corrplot (only for RUSH as originally specified)
  if (name == "RUSH") {
    plot_corrplot(spearman_matrix, "RUSH_2.svg")
  }
}

# ===============================
# Combine All tRNA Datasets
# ===============================
# Ensure these TPM matrices are loaded earlier in your script
# TCGA_tRNA_tpms, GSE110907_tRNA_tpms, etc.

vector_list <- list(
  rownames(TCGA_tRNA_tpms),
  rownames(GSE110907_tRNA_tpms),
  rownames(GSE175462_tRNA_tpms),
  rownames(GSE62182_tRNA_tpms),
  rownames(GSE83527_tRNA_tpms),
  rownames(RUSH_tRNA_tpms)
)
overlap_genes <- Reduce(intersect, vector_list)

all_tRNA <- cbind(
  TCGA_tRNA_tpms[overlap_genes,],
  GSE110907_tRNA_tpms[overlap_genes,],
  GSE175462_tRNA_tpms[overlap_genes,],
  GSE62182_tRNA_tpms[overlap_genes,],
  GSE83527_tRNA_tpms[overlap_genes,],
  RUSH_tRNA_tpms[overlap_genes,]
)

# Combine corresponding sample metadata
all_sample <- rbind(TCGA_sample, GSE110907_sample, GSE175462_sample,
                    GSE62182_sample, GSE83527_sample, RUSH_sample)
write.csv(all_sample, "all_sample.csv")

# ===============================
# Batch Effect Removal (tRNA)
# ===============================
log_tRNA <- log(all_tRNA + 1)
design <- model.matrix(~V2, data = all_sample)
limma_corrected <- removeBatchEffect(log_tRNA, batch = all_sample$batch, design = design)
boxplot(limma_corrected, ylab = "Expression", main = "Batch-Corrected tRNA Data", outline = FALSE)
limma_tRNA_all <- exp(limma_corrected)
stopifnot(all(limma_tRNA_all >= 0))
write.csv(limma_tRNA_all, "limma_tRNA_all.csv")

# ===============================
# Batch Effect Removal (tRF)
# ===============================
vector_list <- list(
  rownames(TCGA_tRF_tpms),
  rownames(GSE110907_tRF_tpms),
  rownames(GSE175462_tRF_tpms),
  rownames(GSE62182_tRF_tpms),
  rownames(GSE83527_tRF_tpms),
  rownames(RUSH_tRF_tpms)
)
overlap_genes <- Reduce(intersect, vector_list)

all_tRF <- cbind(
  TCGA_tRF_tpms[overlap_genes,],
  GSE110907_tRF_tpms[overlap_genes,],
  GSE175462_tRF_tpms[overlap_genes,],
  GSE62182_tRF_tpms[overlap_genes,],
  GSE83527_tRF_tpms[overlap_genes,],
  RUSH_tRF_tpms[overlap_genes,]
)

log_tRF <- log(all_tRF + 1)
design <- model.matrix(~V2, data = all_sample)
limma_corrected_tRF <- removeBatchEffect(log_tRF, batch = all_sample$batch, design = design)
boxplot(limma_corrected_tRF, ylab = "Expression", main = "Batch-Corrected tRF Data", outline = FALSE)
limma_tRF_all <- exp(limma_corrected_tRF)
stopifnot(all(limma_tRF_all >= 0))
write.csv(limma_tRF_all, "limma_tRF_all.csv")

# ===============================
# Correlation Heatmap for Combined Targets
# ===============================
all_filter <- read.csv("all_target_combine.csv", header = TRUE, row.names = 1)
data <- t(all_filter)
cor_matrix_all <- cor(data, method = "spearman")

# Heatmap
plot_correlation_heatmap(cor_matrix_all, "all.svg", "Gene Expression Correlation of All Tissue Samples")

# Corrplot with p-values (optional)
# Uncomment below if needed:
# res_all <- cor.mtest(data, conf.level = 0.95)  # define your own cor.mtest() if needed
# plot_corrplot(cor_matrix_all, "all_corrplot.svg")
