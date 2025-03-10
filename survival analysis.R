#load packages
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(stringr)
library(dbplyr)
library(ggsurvfit)

# getting clinical data for TCGA-LUAD/TCGA-LUSC cohort
clinical_luad<- GDCquery_clinic("TCGA-LUAD")
clinical_lusc<- GDCquery_clinic("TCGA-LUSC")
any(colnames(clinical_luad) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
any(colnames(clinical_lusc) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
which(colnames(clinical_luad) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
which(colnames(clinical_lusc) %in% c("vital_status","days_to_last_follow_up","days_to_death"))
clinical_luad[,c(9,40,46)]
clinical_luad[,c(9,40,45)]

#look at some variables associate with survival
table(clinical_luad$vital_status)
table(clinical_lusc$vital_status)

#change certain values in the way they are encoded
clinical_luad$deceased<-ifelse(clinical_luad$vital_status=="Alive",FALSE, TRUE)
clinical_lusc$deceased<-ifelse(clinical_lusc$vital_status=="Alive",FALSE, TRUE)
table(clinical_luad$deceased)
table(clinical_lusc$deceased)

clinical_luad$overall_survival<- ifelse(clinical_luad$vital_status=="Alive",
                                        clinical_luad$days_to_last_follow_up,
                                        clinical_luad$days_to_death)
clinical_lusc$overall_survival<- ifelse(clinical_lusc$vital_status=="Alive",
                                        clinical_lusc$days_to_last_follow_up,
                                        clinical_lusc$days_to_death)

clinical_luad_extract<- data.frame(clinical_luad$submitter_id,clinical_luad$ajcc_pathologic_stage,clinical_luad$primary_diagnosis,clinical_luad$race,clinical_luad$gender,clinical_luad$ethnicity,clinical_luad$vital_status,clinical_luad$deceased,clinical_luad$overall_survival,clinical_luad$age_at_index,clinical_luad$pack_years_smoked)
clinical_lusc_extract<- data.frame(clinical_lusc$submitter_id,clinical_lusc$ajcc_pathologic_stage,clinical_lusc$primary_diagnosis,clinical_lusc$race,clinical_lusc$gender,clinical_lusc$ethnicity,clinical_lusc$vital_status,clinical_lusc$deceased,clinical_lusc$overall_survival,clinical_lusc$age_at_index,clinical_lusc$pack_years_smoked)
colnames<-c("sample_id","ajcc_pathologic_stage","subtype","race","sex","ethnicity","vital_status","deceased","survival_time","age","smoking_status")
colnames(clinical_luad_extract)<-colnames
colnames(clinical_lusc_extract)<-colnames
clinical_luad_extract$primary_diagnosis<-"LUAD"
clinical_lusc_extract$primary_diagnosis<-"LUSC"

#combine the two meta forms
nsclc_extract<-rbind(clinical_luad_extract,clinical_lusc_extract)
TCGA_rawcounts<-read.csv("TCGA_tRNA.csv",header = T)
TCGA_sample<-read.csv("TCGA_sample.csv",header = T)

to_delete<-rowSums(TCGA_rawcounts==0)/ncol(TCGA_rawcounts)>=0.5
TCGA_rm<-TCGA_rawcounts[!to_delete,] 
row.names(TCGA_sample)<-TCGA_sample$V1
row.names(TCGA_rm)<-TCGA_rm$tRNA_ID
TCGA_rm$tRNA_ID<-NULL
TCGA_sample$V1<-NULL

# Read the TPM normalized TCGA data
TCGA_tpm<-read.csv("TCGA_tRNA_tpms.csv",header = T, row.names = 1)

to_delete<-rowSums(TCGA_tpm==0)/ncol(TCGA_tpm)>=0.5
TCGA_tpm_rm<-TCGA_tpm[!to_delete,]
TCGA_tpm_filter<-TCGA_tpm_rm[rownames(TCGA_tpm_rm) %in% c("tRNA-Val-CAC-2-1", "tRNA-Leu-AAG-2-3", "tRNA-Lys-CTT-3-1", "tRNA-Val-CAC-1-5", "tRNA-Ala-TGC-3-2", "tRNA-Asp-GTC-1-1"), ]

colnames(TCGA_tpm_filter)<- str_extract(colnames(TCGA_tpm_filter), "^[^_]+_[^_]+_[^_]+")
TCGA_extract_match<-nsclc_extract[nsclc_extract$sample_id %in% colnames(TCGA_tpm_filter),]
tpm_filter_transfer<- t(TCGA_tpm_filter)
tpm_filter_transfer<-as.data.frame(tpm_filter_transfer) 
tpm_filter_transfer$sample_id<-rownames(tpm_filter_transfer)
TCGA_metadata_R2_tpm<-merge(tpm_filter_transfer, TCGA_extract_match,by= "sample_id")

#create the Kaplan Meier Curves input csv file for for target gene
#tRNA-Lys-CTT-3-1
svg("tRNA-Lys-CTT-3-1_plot.svg",width = 8, height = 6)

km_curve<-read.csv("tRNA-Lys-CTT-3-1.csv",header = T)
km_curve$compare<-ifelse(km_curve$compare >= 561, 1, 0) #devide the samples based on the cutoff
km_curve$time<-km_curve$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,event)~compare, data = km_curve)
summary(kmcurve)
ggsurvplot(kmcurve, xlim=c(0,72),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = TRUE, surv.median.line = "hv", 
           legend.labs=c("Low tRNA-Lys-CTT-3-1 TPM","High tRNA-Lys-CTT-3-1 TPM"), legend.title="",
           surv.scale = "percent", palette = c("#403990","#cf3d3e"),
           title="",risk.table = TRUE ,risk.table.title="Number at risk")

dev.off()

#tRNA-Val-CAC-2-1
svg("tRNA-Val-CAC-2-1_plot.svg",width = 8, height = 6)
km_curve<-read.csv("tRNA-Val-CAC-2-1.csv",header = T)
km_curve$compare<-ifelse(km_curve$compare >= 370, 1, 0) #devide the samples based on cut-off
km_curve$time<-km_curve$time/30
kmcurve<-survfit(Surv(time,event)~compare, data = km_curve)
ggsurvplot(kmcurve, xlim=c(0,72),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = TRUE,surv.median.line = "hv",
           legend.labs=c("Low tRNA-Val-CAC-2-1 TPM","High tRNA-Val-CAC-2-1 TPM"), legend.title="",
           surv.scale = "percent", palette = c("#403990","#cf3d3e"),
           title="",risk.table = TRUE ,risk.table.title="Number at risk")
dev.off()

#tRNA-Leu-AAG-2-3
svg("tRNA-Leu-AAG-2-3_plot.svg",width = 8, height = 6)
km_curve<-read.csv("tRNA-Leu-AAG-2-3.csv",header = T)
km_curve$compare<-ifelse(km_curve$compare >= 4492, 1, 0) #devide the samples based on cut-off
km_curve$time<-km_curve$time/30
kmcurve<-survfit(Surv(time,event)~compare, data = km_curve)
ggsurvplot(kmcurve, xlim=c(0,72),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = TRUE,surv.median.line = "hv",
           legend.labs=c("Low tRNA-Leu-AAG-2-3 TPM","High tRNA-Leu-AAG-2-3 TPM"), legend.title="",
           surv.scale = "percent", palette = c("#403990","#cf3d3e"),
           title="",risk.table = TRUE ,risk.table.title="Number at risk")

dev.off()

#comparing survival time between groups aka log rank test
logrankres_Leu <- survdiff(Surv(time,event)~compare, data = km_curve)

#assign the risk score for each gene biomarkers
TCGA_metadata_R2_tpm$`tRNA-Leu-AAG-2-3`<-ifelse(TCGA_metadata_R2_tpm$`tRNA-Leu-AAG-2-3`>= 4492, 1, 0)
TCGA_metadata_R2_tpm$`tRNA-Val-CAC-2-1`<-ifelse(TCGA_metadata_R2_tpm$`tRNA-Val-CAC-2-1`>= 370, 1, 0)
TCGA_metadata_R2_tpm$`tRNA-Lys-CTT-3-1`<-ifelse(TCGA_metadata_R2_tpm$`tRNA-Lys-CTT-3-1`>= 561, 1, 0)
TCGA_metadata_R2_tpm$survival_time<-TCGA_metadata_R2_tpm$survival_time/30
TCGA_metadata_R2_tpm$deceased<-ifelse(TCGA_metadata_R2_tpm$deceased == "TRUE", 1, 0)

# Combine the scores into a single signature score
TCGA_metadata_R2_tpm$signature_score <- TCGA_metadata_R2_tpm$`tRNA-Leu-AAG-2-3`+TCGA_metadata_R2_tpm$`tRNA-Val-CAC-2-1`+TCGA_metadata_R2_tpm$`tRNA-Lys-CTT-3-1`
# Find the optimal cut-off for the risk score
opt_cut <- surv_cutpoint(TCGA_metadata_R2_tpm, time = "survival_time", event = "deceased", variables = "signature_score")
opt_cut 

svg("combine_tRNA.svg",width = 8, height = 6)
TCGA_metadata_R2_tpm$risk_group <- ifelse(TCGA_metadata_R2_tpm$signature_score > 1, "High", "Low")

kmcurve<-survfit(Surv(survival_time,deceased)~risk_group, data = TCGA_metadata_R2_tpm)

ggsurvplot(kmcurve, xlim=c(0,72),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = TRUE,surv.median.line = "hv",
           legend.labs=c("High risk score","Low risk score"), legend.title="",
           surv.scale = "percent", palette = c("#fcbb44","#839dd1"),
           title="",risk.table = TRUE ,risk.table.title="Number at risk")

dev.off()

#cox regression

TCGA_metadata_R2_tpm$`smoking_status` <- ifelse(is.na(TCGA_metadata_R2_tpm$`smoking_status`), 0, 1)
TCGA_metadata_R2_tpm<-read.csv("TCGA_meta_tpm.csv",header = T)

TCGA_metadata_R2_tpm$risk_group <- as.factor(TCGA_metadata_R2_tpm$risk_group)
TCGA_metadata_R2_tpm$risk_group <- relevel(TCGA_metadata_R2_tpm$risk_group, ref = "Low")
TCGA_metadata_R2_tpm$sex <- relevel(TCGA_metadata_R2_tpm$sex, ref = "female")
TCGA_metadata_R2_tpm$ajcc_pathologic_stage <- as.factor(TCGA_metadata_R2_tpm$ajcc_pathologic_stage)
TCGA_metadata_R2_tpm$ajcc_pathologic_stage<-relevel(TCGA_metadata_R2_tpm$ajcc_pathologic_stage,ref="Stage I")

cox_model<-coxph(Surv(survival_time,deceased) ~ risk_group +subtype +sex +age +smoking_status+ajcc_pathologic_stage, data =  TCGA_metadata_R2_tpm)
cox_model

library(broom)
library(forestplot)

tidied_cox <- broom::tidy(cox_model)

# Calculate the confidence intervals for the hazard ratios and draw table
tidied_cox$HR <- exp(tidied_cox$estimate)
tidied_cox$lower_CI <- exp(tidied_cox$estimate - 1.96 * tidied_cox$std.error)
tidied_cox$upper_CI <- exp(tidied_cox$estimate + 1.96 * tidied_cox$std.error)


tidied_cox<-read.csv("tidied_cox.csv",header = T)

tidied_cox_2<-cbind(c("Subgroup",tidied_cox$term),
                    c("Number",tidied_cox$number),
                    c("Adjusted HR(95% CI)", tidied_cox$HR),
                    c("P-value",tidied_cox$p.value),
                    c("Significance",tidied_cox$significance))

nd<-rbind("\n",tidied_cox[,c(6:8)])
nd[,c(1:3)]<-lapply(nd[,c(1:3)],as.numeric)
nd$p.value[c(3, 6, 9, 14, 17)] <- 1

svg("forest_plot.svg",width = 10, height = 6)
fig<-forestplot(labeltext = tidied_cox_2,
                nd,
                zero=1,
                graph.pos=4,
                boxsize=0.30,
                lwd.ci=1.7,
                line.margin=unit(1,"mm"),
                lineheight= unit(1,"mm"),
                colgap =unit(5,"mm"),
                xticks= c(0.5,1.0,1.5),
                col=fpColors(box = "#110075",summary = "black",
                             lines = "black",zero = "grey60"),
                txt_gp=fpTxtGp(ticks = gpar(cex=0.7),cex = 0.9),
                hrzl_lines=list("1"=gpar(lty=1),
                                "2"=gpar(lty=1),
                                "21"=gpar(lty=1)),
                graphwidth=unit(0.20,"npc"),
                fn.ci_norm="fpDrawNormalCI",
                is.summary=c(T,T,F,F,T,F,F,T,F,F,T,F,T,F,F,T,F,F,F,F))
plot(fig)

dev.off()
