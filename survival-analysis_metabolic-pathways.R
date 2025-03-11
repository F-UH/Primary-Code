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
nsclc_extract$sample_id<-gsub("-", "_",nsclc_extract$sample_id)
write.csv(nsclc_extract,file = "nsclc_meta.csv")

#read TPM normalized data
TCGA_tpm<- read.table("TCGA_LUSC_LUAD_TPM.txt", sep = "\t", header = TRUE)
rownames(TCGA_tpm)<-TCGA_tpm$ID
TCGA_tpm$ID<-NULL

#remove low counts and filter target genes
to_delete<-rowSums(TCGA_tpm==0)/ncol(TCGA_tpm)>=0.5
TCGA_tpm_rm<-TCGA_tpm[!to_delete,]

#filter the genes of interest
target_gene<- read.csv("Metabolic_pathways.csv", header = T)
gene_list<- target_gene$Symbol
TCGA_target <- TCGA_tpm_rm[rownames(TCGA_tpm_rm) %in% gene_list, ]
TCGA_target <- na.omit(TCGA_target)
colnames(TCGA_target)<- gsub("\\.", "_", colnames(TCGA_target))
colnames(TCGA_target)<- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", colnames(TCGA_target))

TCGA_extract_match<-nsclc_extract[nsclc_extract$sample_id %in% colnames(TCGA_target),]
tpm_target_transfer<- t(TCGA_target)
tpm_target_transfer<-as.data.frame(tpm_target_transfer) 

#calculate Z-scores of the tagret genes
str(tpm_target_transfer)
z_score_data <- apply(tpm_target_transfer, 2, function(column) {
  (column - mean(column, na.rm = TRUE)) / sd(column, na.rm = TRUE)
})
z_score_data<-as.data.frame(z_score_data)
z_score_data$sample_id<-rownames(z_score_data)
TCGA_meta_z<-merge(z_score_data, TCGA_extract_match,by= "sample_id")
write.csv( TCGA_meta_z, "TCGA_meta_z.csv")


#fit the data into a cox model
svg("km_all.svg",width = 8, height = 6)

cox_data<-read.csv("COX_DATA.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
cox_data <- na.omit(cox_data)  
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)
summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  # HR for "Low" group
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#fcbb44","#839dd1"),
           title="Metabolic Pathways",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()

######specfic pathway#####
# Lysine degradation 
#fit the data into a cox model
svg("km_lysine.svg",width = 8, height = 6)
cox_data<-read.csv("lysine_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]


# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)

# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Lysine Degradation",risk.table = TRUE ,risk.table.title="Number at risk")

km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()

# Alanine metabolism
#fit the data into a cox model
svg("km_Alanine.svg",width = 8, height = 6)
cox_data<-read.csv("Alanine_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]

# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Alanine, Aspartate, and Glutamate Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


# Biosynthesis of amino acids
#fit the data into a cox model
svg("km_aa.svg",width = 8, height = 6)
cox_data<-read.csv("AA_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Biosynthesis of Amino Acids",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


# Nucleotide metabolism
#fit the data into a cox model
svg("km_Nucleotide.svg",width = 8, height = 6)
cox_data<-read.csv("Nucleotide_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  # HR for "Low" group
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Nucleotide Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


#Purine metabolism
svg("km_purine.svg",width = 8, height = 6)

cox_data<-read.csv("purine_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  # HR for "Low" group
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Purine Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


#Galactose metabolism
svg("km_galactose.svg",width = 8, height = 6)

cox_data<-read.csv("galactose_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  # HR for "Low" group
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Galactose Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 21, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


#Inositol phosphate metabolism
svg("km_inositol.svg",width = 8, height = 6)

cox_data<-read.csv("inositol_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp") 
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Inositol Phosphate Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()



#Biosynthesis of cofactors
svg("km_cofators.svg",width = 8, height = 6)

cox_data<-read.csv("cofactors_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title="Biosynthesis of Cofactors",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()

#Other types of O-glycan biosynthesis
svg("km_O-glycan.svg",width = 8, height = 6)

cox_data<-read.csv("Oglycan_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  # HR for "Low" group
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title = "Other Types of O-glycan Biosynthesis",risk.table = TRUE ,risk.table.title="Number at risk")

km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


#Taurine and hypotaurine metabolism
svg("km_taurine.svg",width = 8, height = 6)

cox_data<-read.csv("Taurine_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp") 
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3)  
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title = "Taurine and Hypotaurine Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()


#Oxidative phosphorylation
svg("km_oxidative.svg",width = 8, height = 6)

cox_data<-read.csv("Oxidative_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
# Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title = "Oxidative Phosphorylation",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()

#Phosphatidylinositol signaling system
svg("km_Phosphatidylinositol.svg",width = 8, height = 6)

cox_data<-read.csv("Phosphatidylinositol_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp") 
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
#Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title = "Phosphatidylinositol Signaling System",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()

#Carbon metabolism
svg("km_carbon.svg",width = 8, height = 6)

cox_data<-read.csv("Carbon_cox.csv", header = T)
cox_data <- cox_data[, !(colnames(cox_data) %in% "sample_id")]
# Fit Cox model with all genes
cox_data <- na.omit(cox_data)  # Remove rows with any NA values
cox_model <- coxph(Surv(time, status) ~ ., data = cox_data)

summary(cox_model)

#calculate the risk score
cox_data$risk_score <- predict(cox_model, type = "lp")
head(cox_data$risk_score)

# risk stratifycation (Median split)
cox_data$risk_group <- ifelse(cox_data$risk_score > median(cox_data$risk_score), "high", "low")

cox_data$time<-cox_data$time/30 #change the time from days to months

#draw Kaplan Meier curve
kmcurve<-survfit(Surv(time,status)~risk_group, data = cox_data)
summary(kmcurve)
#Cox proportional hazards model
cox_model <- coxph(Surv(time,status)~risk_group, data = cox_data)
cox_summary <- summary(cox_model)
# Extract HR for "Low" risk group (which is the model output)
HR_low <- round(cox_summary$coefficients[1,2], 3) 
HR_high <- round(1 / HR_low, 3)
CI_lower_high <- round(1 / cox_summary$conf.int[1,4], 3)  
CI_upper_high <- round(1 / cox_summary$conf.int[1,3], 3)
cat("HR for High risk group:", HR_high, "(95% CI:", CI_lower_high, "-", CI_upper_high, ")\n")
p_value <- formatC(cox_summary$coefficients[1,5], format = "e", digits = 2)  # Scientific notation for p-value

km_plot<-ggsurvplot(kmcurve, xlim=c(0,120),break.x.by = 12, ylab="Survival probablity",
           xlab="Time", pval = FALSE, surv.median.line = "hv", 
           legend.labs=c("High Risk Group","Low Risk Group"), legend.title="",
           surv.scale = "percent", palette = c("#cf3d3e","#403990"),
           title = "Carbon Metabolism",risk.table = TRUE ,risk.table.title="Number at risk")
km_plot$plot <- km_plot$plot +
  annotate("text", x = 20, y = 0.1, 
           label = paste0( "HR = ", HR_high, " (95% CI: ", CI_lower_high, "–", CI_upper_high, ")\n", "p = ", p_value),
           size = 4,fontface = "plain")
print(km_plot)
dev.off()
