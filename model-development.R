# load packages
library(MetaIntegrator)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(pROC)
library(rstatix)
library(svglite)

getwd()
current_dir <- getwd()
list_files = list.files(current_dir,pattern ="pheno.txt")
file_names <- list_files %>%
  lapply(function(x){strsplit(x,split = "_")[[1]][1]})

  # Load data from the original data source
for (i in seq_along(file_names)) {
  cat("Processing file:", file_names[[i]], "\n")
  cat("Loading original data...\n")
  dataObj1 <- tinyMetaObject$originalData$PBMC.Study.1
  cat("Reading CSV file...\n")
  data <- read.csv(paste0(file_names[[i]], "_design.csv"))
  cat("Reading TXT file and converting to a matrix...\n")
  matrix_data_raw <- read.table(paste0(file_names[[i]], "_batch.txt"))
  matrix_data <- as.matrix(matrix_data_raw)
  cat("Reading phenotype data...\n")
  pheno <- read.table(paste0(file_names[[i]], "_pheno.txt"))
  cat("Converting data to numeric and assigning column names...\n")
  test1 <- as.numeric(data[, 2:dim(data)[2]])
  names(test1) <- names(data)[2:dim(data)[2]]
  cat("Setting 'class' attribute in 'dataObj1'...\n")
  dataObj1$class <- test1
  cat("Setting 'expr' attribute in 'dataObj1' and applying log2 transformation...\n")
  dataObj1$expr <- matrix_data
  dataObj1$expr <- log2(dataObj1$expr + 1)
  cat("Setting 'keys' attribute in 'dataObj1'...\n")
  test1 <- row.names(matrix_data_raw)
  names(test1) <- row.names(matrix_data_raw)
  dataObj1$keys <- test1
  cat("Setting 'pheno' attribute in 'dataObj1'...\n")
  dataObj1$pheno <- pheno
  cat("Setting 'formattedName' attribute in 'dataObj1'...\n")
  dataObj1$formattedName <- file_names[[i]]
  cat("Assigning 'dataObj1' to a variable...\n")
  assign(paste0("dataObj_", i), dataObj1)
  cat("File processed:", file_names[[i]], "\n")
}

file_names

# create metaObjects
discovery_datasets <- list(dataObj_8,dataObj_11,dataObj_12,dataObj_13,dataObj_14)
names(discovery_datasets) = c(dataObj_8$formattedName,dataObj_11$formattedName,dataObj_12$formattedName,dataObj_13$formattedName,dataObj_14$formattedName)
exampleMetaObj=list() 
exampleMetaObj$originalData <- discovery_datasets

validation_datasets <- list(dataObj_15,dataObj_19,dataObj_20,dataObj_21,dataObj_22)
names(validation_datasets) = c(dataObj_15$formattedName,dataObj_19$formattedName, dataObj_20$formattedName,dataObj_21$formattedName,dataObj_22$formattedName)
exampleMetaObj2=list() 
exampleMetaObj2$originalData <- validation_datasets

independent_datasets <- list(dataObj_7)
names(independent_datasets) = c(dataObj_7$formattedName)
exampleMetaObj3<-list()
exampleMetaObj3$originalData <- independent_datasets

checkDataObject(exampleMetaObj, "Meta", "Pre-Analysis")

#warnings()

library(multtest)
library(gtable)

#run metaIntegrator
exampleMetaObj <- runMetaAnalysis(exampleMetaObj, maxCores=1)
exampleMetaObj2 <- runMetaAnalysis(exampleMetaObj2, maxCores=1)
exampleMetaObj3 <- runMetaAnalysis(exampleMetaObj3, maxCores=1)

str(exampleMetaObj, max.level = 2)

##filter results
test1_discovery <- filterGenes(exampleMetaObj, isLeaveOneOut = TRUE, FDRThresh = 0.3, effectSizeThresh = 0.2)
test1_discovery$filterResults$FDR0.3_es0.2_nStudies1_looaTRUE_hetero0
summaryROCPlot(metaObject = exampleMetaObj, filterObject = test1_discovery$filterResults$FDR0.3_es0.2_nStudies1_looaTRUE_hetero0, bootstrapReps = 100)
summary<-summarizeFilterResults(test1_discovery,"FDR0.3_es0.2_nStudies1_looaTRUE_hetero0")
summary_up<-as.data.frame(summary$pos)
summary_down<-as.data.frame(summary$neg)
write.csv(summary_up,"summary_filterresult_pos.csv")
write.csv(summary_down,"summary_filterresult_neg.csv")

#Hold-out validation
summaryROCPlot(metaObject = exampleMetaObj2, filterObject = test1_discovery$filterResults$FDR0.3_es0.2_nStudies1_looaTRUE_hetero0, bootstrapReps = 100)

#Independent validation by RUSH
rocPlot(test1_discovery$filterResults$FDR0.3_es0.2_nStudies1_looaTRUE_hetero0, dataObj_7, title = "ROC plot for independent validation")

# forwardRes
forwardRes <- forwardSearch(test1_discovery,
                            test1_discovery$filterResults$FDR0.3_es0.2_nStudies1_looaTRUE_hetero0,
                            forwardThresh = 0)

forwardRes

# ROC
svg("Discorvey_ROC.svg", width=7, height=5) 
summary<-summaryROCPlot(metaObject = exampleMetaObj, filterObject = forwardRes, bootstrapReps = 100)
dev.off()


# Heatmap of effect across studies
svg("heatmap_trna.svg", width=6, height=4)
heatmapPlot(exampleMetaObj_all, forwardRes, colorRange = c(-1, 1), geneOrder = FALSE, datasetOrder = FALSE, displayPooled = F,useFormattedNames = TRUE)
dev.off()

# forest plot for each gene
svg("forest_tRNA-Val-CAC-2-1.svg", width=5, height=6) 
forestPlot(exampleMetaObj_all, geneName="tRNA-Val-CAC-2-1",  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#e64b35ff", textColor = "#3c5488ff")
dev.off()
svg("forest_tRNA-Leu-AAG-2-3.svg", width=5, height=6)
forestPlot(exampleMetaObj_all, geneName="tRNA-Leu-AAG-2-3",  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#e64b35ff", textColor = "#3c5488ff")
dev.off()
svg("forest_tRNA-Lys-CTT-3-1.svg", width=5, height=6)
forestPlot(exampleMetaObj_all, geneName="tRNA-Lys-CTT-3-1",  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#e64b35ff", textColor = "#3c5488ff")
dev.off()
svg("forest_tRNA-Val-CAC-1-5.svg", width=5, height=6)
forestPlot(exampleMetaObj_all, geneName="tRNA-Val-CAC-1-5",  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#e64b35ff", textColor = "#3c5488ff")
dev.off()
svg("forest_tRNA-Ala-TGC-3-2.svg", width=5, height=6)
forestPlot(exampleMetaObj_all, geneName="tRNA-Ala-TGC-3-2",  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#5d9ad3", textColor = "#3c5488ff")
dev.off()
svg("forest_tRNA-Asp-GTC-1-1.svg", width=5, height=6)
forestPlot(exampleMetaObj_all, geneName="tRNA-Asp-GTC-1-1" ,  boxColor = "#00A087FF", whiskerColor = "darkgrey", 
           zeroLineColor = "black", summaryColor = "#5d9ad3", textColor = "#3c5488ff")
dev.off()

# ROC of Hold-out validation
svg("Holdout_ROC.svg", width=7, height=5) 
summaryROCPlot(metaObject = exampleMetaObj2, filterObject = forwardRes, bootstrapReps = 100)
dev.off()

# ROC of Independent validation
svg("Independet_ROC.svg", width=7, height=5) 
rocPlotValidation <- rocPlot(forwardRes, dataObj_7, title = "Independent validation in RUSH cohort")
dev.off()


# draw a box plot based on the T-score of discovery phase (training)
GSE110907_score<-as.data.frame(calculateScore(forwardRes, dataObj_8))
colnames(GSE110907_score)<- "Tscore"
GSE110907_pheno <- dataObj_8$pheno
GSE110907_score$group <- GSE110907_pheno$group
GSE110907_score$dataset <- "GSE110907"

GSE175462_score<-as.data.frame(calculateScore(forwardRes, dataObj_11))
colnames(GSE175462_score)<- "Tscore"
GSE175462_pheno <- dataObj_11$pheno
GSE175462_score$group <- GSE175462_pheno$group
GSE175462_score$dataset <- "GSE175462"

GSE62182_score<-as.data.frame(calculateScore(forwardRes, dataObj_12))
colnames(GSE62182_score)<- "Tscore"
GSE62182_pheno <- dataObj_12$pheno
GSE62182_score$group <- GSE62182_pheno$group
GSE62182_score$dataset <- "GSE62182"

GSE83527_score<-as.data.frame(calculateScore(forwardRes, dataObj_13))
colnames(GSE83527_score)<- "Tscore"
GSE83527_pheno <- dataObj_13$pheno
GSE83527_score$group <- GSE83527_pheno$group
GSE83527_score$dataset <- "GSE83527"

TCGA_score<-as.data.frame(calculateScore(forwardRes, dataObj_14))
colnames(TCGA_score)<- "Tscore"
TCGA_pheno <- dataObj_14$pheno
TCGA_score$group <- TCGA_pheno$group
TCGA_score$dataset <- "TCGA"

bind_t<-bind_rows(GSE110907_score,GSE175462_score, GSE62182_score, GSE83527_score, TCGA_score)
bind_t$diagnosis <- ifelse(bind_t$group == "tumor", "1", "0")


# Calculate p value by Mann-Whitney U test
pval <- bind_t %>%
  group_by(dataset) %>%
  wilcox_test(Tscore ~ group) %>% 
  add_significance(p.col = 'p', cutpoints = c(0,0.001,0.01,0.05,1), symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='dataset')

# Draw box plot by ggplot
box_plot <- ggplot() +
  geom_boxplot(data = bind_t, mapping = aes(x = dataset, y = Tscore, fill = group),
               width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.65)) +
  scale_fill_manual(values = c('#80c5a2', '#f27873')) +
  stat_pvalue_manual(pval, label = 'p.signif', tip.length = 0.03) +
  labs(x = 'Discovery Datasets', y = 'T-score') +
  guides(fill = guide_legend(title = 'Diagnosis')) +
  theme(
    axis.text = element_text(color = 'black'),
    legend.position = c(0.76, 0.1),
    legend.direction = 'horizontal',
    legend.background = element_rect(fill = "white", colour = "white"), 
    legend.key = element_rect(fill = "white", colour = "white"), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    panel.background = element_rect(fill = "white", colour = "grey"), 
    panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"), 
    axis.line = element_line(colour = "black") 
  ) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1), expand = c(0, 0))

box_plot

ggsave("Discovery_boxplot_1.svg", plot = box_plot, width=6, height=3)

#draw a box plot based on the T-score of hold-out validation 
vGSE110907_score<-as.data.frame(calculateScore(forwardRes, dataObj_15))
colnames(vGSE110907_score)<- "Tscore"
vGSE110907_pheno <- dataObj_15$pheno
vGSE110907_score$group <- vGSE110907_pheno$group
vGSE110907_score$dataset <- "GSE110907"

vGSE175462_score<-as.data.frame(calculateScore(forwardRes, dataObj_19))
colnames(vGSE175462_score)<- "Tscore"
vGSE175462_pheno <- dataObj_19$pheno
vGSE175462_score$group <- vGSE175462_pheno$group
vGSE175462_score$dataset <- "GSE175462"

vGSE62182_score<-as.data.frame(calculateScore(forwardRes, dataObj_20))
colnames(vGSE62182_score)<- "Tscore"
vGSE62182_pheno <- dataObj_20$pheno
vGSE62182_score$group <- vGSE62182_pheno$group
vGSE62182_score$dataset <- "GSE62182"

vGSE83527_score<-as.data.frame(calculateScore(forwardRes, dataObj_21))
colnames(vGSE83527_score)<- "Tscore"
vGSE83527_pheno <- dataObj_21$pheno
vGSE83527_score$group <- vGSE83527_pheno$group
vGSE83527_score$dataset <- "GSE83527"

vTCGA_score<-as.data.frame(calculateScore(forwardRes, dataObj_22))
colnames(vTCGA_score)<- "Tscore"
vTCGA_pheno <- dataObj_22$pheno
vTCGA_score$group <- vTCGA_pheno$group
vTCGA_score$dataset <- "TCGA"

bind<-bind_rows(vGSE110907_score,vGSE175462_score, vGSE62182_score, vGSE83527_score, vTCGA_score)


pval <- bind %>%
  group_by(dataset) %>%
  wilcox_test(Tscore ~ group) %>% 
  add_significance(p.col = 'p', cutpoints = c(0,0.001,0.01,0.05,1), symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='dataset')

box_plot <- ggplot() +
  geom_boxplot(data = bind, mapping = aes(x = dataset, y = Tscore, fill = group),
               width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.65)) +
  scale_fill_manual(values = c('#80c5a2', '#f27873')) +
  stat_pvalue_manual(pval, label = 'p.signif', tip.length = 0.03) +
  labs(x = 'Hold-out Validation Datasets', y = 'T-score') +
  guides(fill = guide_legend(title = 'Diagnosis')) +
  theme(
    axis.text = element_text(color = 'black'),
    legend.position = c(0.76, 0.1),
    legend.direction = 'horizontal',
    legend.background = element_rect(fill = "white", colour = "white"), 
    legend.key = element_rect(fill = "white", colour = "white"),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    panel.background = element_rect(fill = "white", colour = "grey"), 
    panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"), 
    axis.line = element_line(colour = "black") 
  ) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1), expand = c(0, 0))

box_plot 

ggsave("Hold-out Validation_boxplot_1.svg", plot = box_plot, width=6, height=3)

#draw a box plot based on the T-score of independent validation 
RUSH_score<-as.data.frame(calculateScore(forwardRes, dataObj_7))
colnames(RUSH_score)<- "Tscore"
RUSH_pheno <- dataObj_7$pheno
RUSH_score$group <- RUSH_pheno$group
RUSH_score$dataset <- "RUSH"
write.csv(RUSH_score,file = "RUSH_score.csv")

pval <- RUSH_score %>%
  group_by(dataset) %>%
  wilcox_test(Tscore ~ group) %>% 
  add_significance(p.col = 'p', cutpoints = c(0,0.001,0.01,0.05,1), symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='dataset')


box_plot <- ggplot() +
  geom_boxplot(data = RUSH_score, mapping = aes(x = dataset, y = Tscore, fill = group),
               width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.65)) +
  scale_fill_manual(values = c('#80c5a2', '#f27873')) +
  stat_pvalue_manual(pval, label = 'p.signif', tip.length = 0.03) +
  labs(x = 'Independent Validation Datasets', y = 'T-score') +
  guides(fill = guide_legend(title = 'Diagnosis')) +
  theme(
    axis.text = element_text(color = 'black'),
    legend.position = c(0.7, 0.2),
    legend.direction = 'vertical',
    legend.background = element_rect(fill = "white", colour = "white"), 
    legend.key = element_rect(fill = "white", colour = "white"), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(), 
    panel.background = element_rect(fill = "white", colour = "grey"), 
    panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"), 
    axis.line = element_line(colour = "black") 
  ) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1), expand = c(0, 0))

box_plot 

ggsave("Independent Validation_boxplot_1.svg", plot = box_plot, width=3, height=3)

