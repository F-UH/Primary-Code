# Load packages
library(dplyr)
library(stringr)
library(DESeq2)
library(limma)
library(ggplot2)

# Read files
TCGA_tpm<- read.csv("TCGA_tpm.csv",header = 1, row.names = 1)
TCGA_sample<-read.csv("TCGA_sample.csv",header = F)
TCGA_sample$batch<-"TCGA"
TCGA_sample$type<-"tissue"
GSE62182_tpm<- read.csv("GSE62182_tpm.csv",header = 1, row.names = 1)
GSE62182_sample<-read.csv("GSE62182_filter_sample.csv",header = F)
GSE62182_sample$batch<-"GSE62182"
GSE62182_sample$type<-"tissue"
GSE83527_tpm<- read.csv("GSE83527_tpm.csv",header = 1, row.names = 1)
GSE83527_sample<-read.csv("GSE83527_sample.csv",header = F)
GSE83527_sample$batch<-"GSE83527"
GSE83527_sample$type<-"tissue"
GSE175462_tpm<- read.csv("GSE175462_tpm.csv",header = 1, row.names = 1)
GSE175462_sample<-read.csv("GSE175462_sample.csv",header = F)
GSE175462_sample$batch<-"GSE175462"
GSE175462_sample$type<-"tissue"
GSE110907_tpm<- read.csv("GSE110907_tpm.csv",header = 1, row.names = 1)
GSE110907_sample<-read.csv("GSE110907_sample.csv",header = F)
GSE110907_sample$batch<-"GSE110907"
GSE110907_sample$type<-"tissue"
RUSH_tpm<-read.csv("RUSH_tpm.csv",header = T, row.names = 1)
RUSH_sample<-read.csv("RUSH_sample.csv",header = F)
RUSH_sample$batch<-"RUSH"
RUSH_sample$type<-"exosome"

#remove low counts
to_delete<-rowSums(TCGA_tpm==0)/ncol(TCGA_tpm)>=0.5
TCGA_rm<-TCGA_tpm[!to_delete,]
to_delete<-rowSums(GSE62182_tpm==0)/ncol(GSE62182_tpm)>=0.5
GSE62182_rm<-GSE62182_tpm[!to_delete,]
to_delete<-rowSums(GSE83527_tpm==0)/ncol(GSE83527_tpm)>=0.5
GSE83527_rm<-GSE83527_tpm[!to_delete,]
to_delete<-rowSums(GSE175462_tpm==0)/ncol(GSE175462_tpm)>=0.5
GSE175462_rm<-GSE175462_tpm[!to_delete,]
to_delete<-rowSums(GSE110907_tpm==0)/ncol(GSE110907_tpm)>=0.5
GSE110907_rm<-GSE110907_tpm[!to_delete,]
to_delete<-rowSums(RUSH_tpm==0)/ncol(RUSH_tpm)>=0.5
RUSH_rm<-RUSH_tpm[!to_delete,]

#combine data
all_tRNA<-cbind(GSE83527_rm[overlap_PCA,],GSE62182_rm[overlap_PCA,],GSE175462_rm[overlap_PCA,],GSE110907_rm[overlap_PCA,],TCGA_rm[overlap_PCA,],RUSH_rm[overlap_PCA,]) 
all_sample<-rbind(GSE83527_sample,GSE62182_sample,GSE175462_sample,GSE110907_sample,TCGA_sample,RUSH_sample)   
all_tRNA_matrix<-t(all_tRNA)

#remove_batch effect & PCA
design<-model.matrix(~V2,data=all_sample) 
limma<-removeBatchEffect(all_tRNA,batch=all_sample$batch,design = design) 
limma.data<-t(limma)
limma_matrix<-apply(limma.data, 2,as.numeric)
any(is.na(limma_matrix)) 
scaled_limma<-scale(limma_matrix)
pc<-prcomp(scaled_limma,center=T,scale=T)
print(pc)
plot(pc)
pc$x[,1:2]
pc_results <- data.frame(pc$x[, 1:2], stringsAsFactors = FALSE, group = all_sample$V2, batch = all_sample$batch)

# Plot PCA result
pc_results$batch <- factor(pc_results$batch, levels = c("TCGA", "GSE83527", "GSE62182", "GSE175462", "GSE110907", "RUSH"))
svg("PCA_all_limma.svg",width = 7, height = 5)
pca_plot <- ggplot(pc_results, aes(x = PC1, y = PC2, color = group,shape=factor(all_sample$batch))) +
  geom_point(size = 1.5) +
  stat_ellipse(aes(group = group), geom = "polygon", alpha = 0.2, fill = NA) +
  xlab("PC 1") +
  ylab("PC 2") +
  labs(title = "PCA Analysis of tRNA Datasets",
       shape = "Batch",
       color = "Diagnosis") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        legend.title = element_text(face="bold"),
        plot.title = element_text(face="bold", size=18),
        axis.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=14),
        axis.line = element_line(color="black"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid"),
        legend.key = element_blank(),
        panel.grid.major = element_line(color="grey80"),
        panel.grid.minor = element_blank()) +
  scale_color_brewer(palette="Set2", guide=guide_legend(title="Group")) +
  scale_shape_manual(values=c(19, 15, 17, 18, 20, 13), guide=guide_legend(title="Batch")) +
  guides(color=guide_legend(override.aes=list(size=3))) +
  coord_cartesian(xlim = c(-20, 40))
pca_plot
dev.off()
