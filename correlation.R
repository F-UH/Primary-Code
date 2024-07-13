#Correlation
library(pheatmap)
library(RColorBrewer)

RUSH_filter<-read.csv(file = "RUSH_filter.csv", header = T, row.names = 1)
data<-t(RUSH_filter)
spearman_RUSH<-cor(data,method = "spearman")
pheatmap(spearman_RUSH, frontsize_row = 18, main = "Gene Expression Correlation of RUSH")
colors <- colorRampPalette(brewer.pal(15, "RdBu"))(255)

svg("RUSH.svg", width=10, height=8) 
pheatmap(spearman_RUSH, 
         col = colors,
         fontsize_row = 10,  
         fontsize_col = 10,
         main = "Gene Expression Correlation of RUSH",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)
dev.off() 

library(corrplot)
svg("RUSH_2.svg",width = 10,height = 8)
colors <- colorRampPalette(brewer.pal(10, "RdYlBu"))(255)
all<- corrplot(spearman_RUSH,
               method = "shade",
               type = "upper",
               addCoef.col = "white",
               number.cex = 0.7,
               col = colors,
               tl.cex = 0.8,
               tl.col = "black",
               tl.srt = 45)
dev.off()

TCGA_filter<-read.csv(file = "TCGA_filter.csv", header = T, row.names = 1)
data<-t(TCGA_filter)
spearman_TCGA<-cor(data,method = "spearman")
pheatmap(spearman_TCGA, frontsize_row = 18, main = "Gene Expression Correlation of TCGA")
colors <- colorRampPalette(brewer.pal(15, "RdBu"))(255)
svg("TCGA.svg", width=10, height=8)
pheatmap(spearman_TCGA, 
         col = colors,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Gene Expression Correlation of TCGA",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE) 
dev.off()

GSE110907_filter<-read.csv(file = "GSE110907_filter.csv", header = T, row.names = 1)
data<-t(GSE110907_filter)
spearman_GSE110907<-cor(data,method = "spearman")
colors <- colorRampPalette(brewer.pal(15, "RdBu"))(255)
svg("GSE110907.svg",width = 10, height = 8)
pheatmap(spearman_GSE110907, 
         col = colors,
         fontsize_row = 10, 
         fontsize_col = 10,
         main = "Gene Expression Correlation of GSE110907",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)
dev.off()


GSE62182_filter<-read.csv(file = "GSE62182_filter.csv", header = T, row.names = 1)
data<-t(GSE62182_filter)
spearman_GSE62182<-cor(data,method = "spearman")
pheatmap(spearman_GSE62182, frontsize_row = 18, main = "Gene Expression Correlation of GSE62182")
svg("GSE62182.svg",width = 10, height = 8)
pheatmap(spearman_GSE62182, 
         col = colors,
         fontsize_row = 10,  
         fontsize_col = 10,
         main = "Gene Expression Correlation of GSE62182",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)
dev.off()

GSE83527_filter<-read.csv(file = "GSE83527_filter.csv", header = T, row.names = 1)
data<-t(GSE83527_filter)
spearman_GSE83527<-cor(data,method = "spearman")
pheatmap(spearman_GSE83527, frontsize_row = 18, main = "Gene Expression Correlation of GSE83527")
svg("GSE83527.svg", width = 10, height = 8)
pheatmap(spearman_GSE83527, 
         col = colors,
         fontsize_row = 10, 
         fontsize_col = 10,
         main = "Gene Expression Correlation of GSE83527",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)
dev.off()

GSE175462_filter<-read.csv(file = "GSE175462_filter.csv", header = T, row.names = 1)
data<-t(GSE175462_filter)
spearman_GSE175462<-cor(data,method = "spearman")
pheatmap(spearman_GSE175462, frontsize_row = 18, main = "Gene Expression Correlation of GSE175462")
svg("GSE175462.svg",width = 10,height = 8)
pheatmap(spearman_GSE175462, 
         col = colors,
         fontsize_row = 10,  
         fontsize_col = 10,
         main = "Gene Expression Correlation of GSE175462",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)

dev.off()


# Combine all the tissue sample together
RUSH_sample$batch <- "RUSH"
TCGA_sample$batch <- "TCGA"
GSE110907_sample$batch <- "GSE110907"
GSE175462_sample$batch <- "GSE175462"
GSE62182_sample$batch <- "GSE62182"
GSE83527_sample$batch <- "GSE83527"

vector_list <- list(
  rownames(TCGA_tRNA_tpms),
  rownames(GSE110907_tRNA_tpms),
  rownames(GSE175462_tRNA_tpms),
  rownames(GSE62182_tRNA_tpms),
  rownames(GSE83527_tRNA_tpms)
)
overlap_PCA <- Reduce(intersect, vector_list)
all_tRNA<-cbind(TCGA_tRNA_tpms[overlap_PCA,],GSE110907_tRNA_tpms[overlap_PCA,],GSE175462_tRNA_tpms[overlap_PCA,],GSE62182_tRNA_tpms[overlap_PCA,],GSE83527_tRNA_tpms[overlap_PCA,])
all_sample<-rbind(TCGA_sample,GSE110907_sample,GSE175462_sample,GSE62182_sample,GSE83527_sample)
write.csv(all_sample,file = "all_sample.csv")
log_all <- log(all_tRNA+1)
boxplot(log_all,ylab="Data", main="With out batch effect correction", outline=FALSE, notch=FALSE)

#remove batch effect
library(limma)
design<-model.matrix(~V2, data = all_sample)
limma<-removeBatchEffect(log_all,batch = all_sample$batch, design = design)
boxplot(limma,ylab="Data", main="Data after batch effect correction", outline=FALSE, notch=FALSE)

# Applying the inverse log transformation to revert back to counts
limma_tRNA_all <- exp(limma)
any(limma_tRNA_all < 0)
write.csv(limma_tRNA_all, file = "limma_tRNA_all.csv")

#for tRF
vector_list <- list(
  rownames(TCGA_tRF_tpms),
  rownames(GSE110907_tRF_tpms),
  rownames(GSE175462_tRF_tpms),
  rownames(GSE62182_tRF_tpms),
  rownames(GSE83527_tRF_tpms)
)
overlap_PCA <- Reduce(intersect, vector_list)
all_tRF<-cbind(TCGA_tRF_tpms[overlap_PCA,],GSE110907_tRF_tpms[overlap_PCA,],GSE175462_tRF_tpms[overlap_PCA,],GSE62182_tRF_tpms[overlap_PCA,],GSE83527_tRF_tpms[overlap_PCA,])
all_sample<-rbind(TCGA_sample,GSE110907_sample,GSE175462_sample,GSE62182_sample,GSE83527_sample)
log_all <- log(all_tRF+1)
boxplot(log_all,ylab="Data", main="With out batch effect correction", outline=FALSE, notch=FALSE)

design<-model.matrix(~V2, data = all_sample)
limma<-removeBatchEffect(log_all,batch = all_sample$batch, design = design)
boxplot(limma,ylab="Data", main="Data after batch effect correction", outline=FALSE, notch=FALSE)

limma_tRF_all <- exp(limma)
any(limma_tRNA_all < 0)
write.csv(limma_tRF_all, file = "limma_tRF_all.csv")

#draw heatmap
all_filter<-read.csv(file = "all_target_combine.csv", header = T, row.names = 1)
data<-t(all_filter)
cor_mat_all<-cor(data)
pheatmap(cor_mat_all, frontsize_row = 18, main = "Gene Expression Correlation of tissue samples")
colors <- colorRampPalette(brewer.pal(20, "RdBu"))(255)
svg("all.svg",width = 10,height = 8)
pheatmap(spearman_all, 
         col = colors,
         fontsize_row = 10, 
         fontsize_col = 10,
         main = "Gene Expression Correlation of all tissue samples",
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         legend = TRUE,
         breaks = seq(-1, 1, length.out = length(colors)),
         show_rownames = TRUE, show_colnames = TRUE)
dev.off()
