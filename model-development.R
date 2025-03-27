# =====================================
# Load Required Packages
# =====================================
library(MetaIntegrator)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(pROC)
library(rstatix)
library(svglite)
library(multtest)
library(gtable)

# =====================================
# Load and Prepare MetaIntegrator Datasets
# =====================================
current_dir <- getwd()
file_names <- list.files(current_dir, pattern = "pheno.txt") %>%
  strsplit("_") %>%
  sapply(`[`, 1)

for (i in seq_along(file_names)) {
  cat("Processing:", file_names[[i]], "\n")
  dataObj <- tinyMetaObject$originalData$PBMC.Study.1
  
  expr_raw <- read.table(paste0(file_names[[i]], "_batch.txt"))
  dataObj$expr <- log2(as.matrix(expr_raw) + 1)
  
  design <- read.csv(paste0(file_names[[i]], "_design.csv"))
  dataObj$class <- as.numeric(design[, 2:ncol(design)])
  names(dataObj$class) <- names(design)[2:ncol(design)]
  
  dataObj$pheno <- read.table(paste0(file_names[[i]], "_pheno.txt"))
  rownames <- rownames(expr_raw)
  dataObj$keys <- setNames(rownames, rownames)
  
  dataObj$formattedName <- file_names[[i]]
  assign(paste0("dataObj_", i), dataObj)
}

# =====================================
# Create MetaObjects
# =====================================
# Discovery Datasets
discovery_datasets <- list(dataObj_8, dataObj_11, dataObj_12, dataObj_13, dataObj_14)
names(discovery_datasets) <- sapply(discovery_datasets, `[[`, "formattedName")
exampleMetaObj <- list(originalData = discovery_datasets)

# Validation Datasets
validation_datasets <- list(dataObj_15, dataObj_19, dataObj_20, dataObj_21, dataObj_22)
names(validation_datasets) <- sapply(validation_datasets, `[[`, "formattedName")
exampleMetaObj2 <- list(originalData = validation_datasets)

# Independent Dataset
independent_datasets <- list(dataObj_7)
names(independent_datasets) <- sapply(independent_datasets, `[[`, "formattedName")
exampleMetaObj3 <- list(originalData = independent_datasets)

# =====================================
# Meta-Analysis
# =====================================
checkDataObject(exampleMetaObj, "Meta", "Pre-Analysis")

exampleMetaObj <- runMetaAnalysis(exampleMetaObj, maxCores = 1)
exampleMetaObj2 <- runMetaAnalysis(exampleMetaObj2, maxCores = 1)
exampleMetaObj3 <- runMetaAnalysis(exampleMetaObj3, maxCores = 1)

# =====================================
# Filter Genes from Discovery
# =====================================
filterObj <- filterGenes(exampleMetaObj, isLeaveOneOut = TRUE, FDRThresh = 0.3, effectSizeThresh = 0.2)
filtered <- filterObj$filterResults[[1]]

summary <- summarizeFilterResults(filterObj, names(filterObj$filterResults)[1])
write.csv(as.data.frame(summary$pos), "summary_filterresult_pos.csv")
write.csv(as.data.frame(summary$neg), "summary_filterresult_neg.csv")

# =====================================
# ROC & Heatmaps
# =====================================
svg("Discovery_ROC.svg", width = 7, height = 5)
summaryROCPlot(exampleMetaObj, filterObject = filtered, bootstrapReps = 100)
dev.off()

summaryROCPlot(exampleMetaObj2, filterObject = filtered, bootstrapReps = 100)
rocPlot(filtered, dataObj_7, title = "Independent ROC")

# Forward Search Refinement
forwardRes <- forwardSearch(filterObj, filtered, forwardThresh = 0)

svg("Discovery_ROC_forward.svg", width = 7, height = 5)
summaryROCPlot(exampleMetaObj, filterObject = forwardRes, bootstrapReps = 100)
dev.off()

svg("heatmap_trna.svg", width = 6, height = 4)
heatmapPlot(exampleMetaObj, forwardRes, colorRange = c(-1, 1), displayPooled = FALSE, useFormattedNames = TRUE)
dev.off()

# =====================================
# Forest Plots
# =====================================
genes <- c("tRNA-Val-CAC-2-1", "tRNA-Leu-AAG-2-3", "tRNA-Lys-CTT-3-1",
           "tRNA-Val-CAC-1-5", "tRNA-Ala-TGC-3-2", "tRNA-Asp-GTC-1-1")

for (gene in genes) {
  svg(paste0("forest_", gene, ".svg"), width = 5, height = 6)
  forestPlot(exampleMetaObj, geneName = gene,
             boxColor = "#00A087FF", whiskerColor = "darkgrey",
             zeroLineColor = "black", summaryColor = "#e64b35ff", textColor = "#3c5488ff")
  dev.off()
}

# =====================================
# ROC for Validation Sets
# =====================================
svg("Holdout_ROC.svg", width = 7, height = 5)
summaryROCPlot(exampleMetaObj2, filterObject = forwardRes, bootstrapReps = 100)
dev.off()

svg("Independent_ROC.svg", width = 7, height = 5)
rocPlot(forwardRes, dataObj_7, title = "Independent Validation in RUSH")
dev.off()

# =====================================
# Boxplot Function for T-scores
# =====================================
plot_box_scores <- function(datasets, labels, filename, title) {
  scores <- bind_rows(lapply(seq_along(datasets), function(i) {
    score_df <- as.data.frame(calculateScore(forwardRes, datasets[[i]]))
    colnames(score_df) <- "Tscore"
    score_df$group <- datasets[[i]]$pheno$group
    score_df$dataset <- labels[i]
    return(score_df)
  }))
  
  pval <- scores %>%
    group_by(dataset) %>%
    wilcox_test(Tscore ~ group) %>%
    add_significance(p.col = 'p', cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns')) %>%
    add_xy_position(x = 'dataset')
  
  plot <- ggplot(scores, aes(x = dataset, y = Tscore, fill = group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.65)) +
    scale_fill_manual(values = c('#80c5a2', '#f27873')) +
    stat_pvalue_manual(pval, label = 'p.signif', tip.length = 0.03) +
    labs(x = title, y = 'T-score') +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text = element_text(color = 'black'),
          panel.border = element_rect(colour = "black", fill = NA)) +
    scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))
  
  ggsave(filename, plot, width = 6, height = 3)
  return(plot)
}

# Discovery
discovery_labels <- c("GSE110907", "GSE175462", "GSE62182", "GSE83527", "TCGA")
discovery_objs <- list(dataObj_8, dataObj_11, dataObj_12, dataObj_13, dataObj_14)
plot_box_scores(discovery_objs, discovery_labels, "Discovery_boxplot_1.svg", "Discovery Datasets")

# Hold-out
validation_labels <- c("GSE110907", "GSE175462", "GSE62182", "GSE83527", "TCGA")
validation_objs <- list(dataObj_15, dataObj_19, dataObj_20, dataObj_21, dataObj_22)
plot_box_scores(validation_objs, validation_labels, "Holdout_boxplot_1.svg", "Hold-out Validation")

# Independent
indep_score <- as.data.frame(calculateScore(forwardRes, dataObj_7))
colnames(indep_score) <- "Tscore"
indep_score$group <- dataObj_7$pheno$group
indep_score$dataset <- "RUSH"
write.csv(indep_score, "RUSH_score.csv")

indep_pval <- indep_score %>%
  group_by(dataset) %>%
  wilcox_test(Tscore ~ group) %>%
  add_significance(p.col = 'p', cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns')) %>%
  add_xy_position(x = 'dataset')

indep_plot <- ggplot(indep_score, aes(x = dataset, y = Tscore, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.65)) +
  scale_fill_manual(values = c('#80c5a2', '#f27873')) +
  stat_pvalue_manual(indep_pval, label = 'p.signif', tip.length = 0.03) +
  labs(x = 'Independent Validation', y = 'T-score') +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text = element_text(color = 'black'),
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))

ggsave("Independent_boxplot_1.svg", indep_plot, width = 3, height = 3)
