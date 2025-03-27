#### 00. Setup ####
# Load required packages
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(stringr)
library(ggsurvfit)

#### 01. Clinical Data Processing ####
# Download clinical data
clinical_luad <- GDCquery_clinic("TCGA-LUAD")
clinical_lusc <- GDCquery_clinic("TCGA-LUSC")

# Create consistent survival variables
clinical_luad$deceased <- clinical_luad$vital_status != "Alive"
clinical_lusc$deceased <- clinical_lusc$vital_status != "Alive"

clinical_luad$overall_survival <- ifelse(clinical_luad$deceased,
                                         clinical_luad$days_to_death,
                                         clinical_luad$days_to_last_follow_up)
clinical_lusc$overall_survival <- ifelse(clinical_lusc$deceased,
                                         clinical_lusc$days_to_death,
                                         clinical_lusc$days_to_last_follow_up)

# Extract and unify clinical features
cols <- c("submitter_id", "ajcc_pathologic_stage", "primary_diagnosis",
          "race", "gender", "ethnicity", "vital_status", "deceased",
          "overall_survival", "age_at_index", "pack_years_smoked")
clinical_luad_extract <- clinical_luad[, cols]
clinical_lusc_extract <- clinical_lusc[, cols]

colnames(clinical_luad_extract) <- colnames(clinical_lusc_extract) <- c("sample_id", "ajcc_stage", "subtype", 
    "race", "sex", "ethnicity", "vital_status", "deceased", 
    "survival_time", "age", "smoking_status")

clinical_luad_extract$primary_diagnosis <- "LUAD"
clinical_lusc_extract$primary_diagnosis <- "LUSC"

# Merge both cancer types
nsclc_clinical <- rbind(clinical_luad_extract, clinical_lusc_extract)
nsclc_clinical$sample_id <- gsub("-", "_", nsclc_clinical$sample_id)
write.csv(nsclc_clinical, "nsclc_meta.csv", row.names = FALSE)

#### 02. TPM Expression Data and Target Genes ####
# Read expression matrix
tpm <- read.table("TCGA_LUSC_LUAD_TPM.txt", sep = "\t", header = TRUE, row.names = 1)

# Remove low expression genes
tpm_filtered <- tpm[rowSums(tpm == 0) / ncol(tpm) < 0.5, ]

# Filter for target genes
target_gene <- read.csv("Metabolic_pathways.csv")
tpm_target <- tpm_filtered[rownames(tpm_filtered) %in% target_gene$Symbol, ]
tpm_target <- na.omit(tpm_target)

# Clean sample names
colnames(tpm_target) <- gsub("\\.", "_", colnames(tpm_target))
colnames(tpm_target) <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", colnames(tpm_target))

# Match with clinical
matched_clinical <- nsclc_clinical[nsclc_clinical$sample_id %in% colnames(tpm_target), ]
tpm_df <- as.data.frame(t(tpm_target))
tpm_df$sample_id <- rownames(tpm_df)

# Calculate Z-scores
z_data <- as.data.frame(scale(tpm_df[, -ncol(tpm_df)]))
z_data$sample_id <- tpm_df$sample_id

# Merge with matched clinical data
merged_data <- merge(z_data, matched_clinical, by = "sample_id")
write.csv(merged_data, "TCGA_meta_z.csv", row.names = FALSE)

#### 03. Cox Model + KM Plot Function ####
plot_km_survival <- function(input_csv, output_svg, title, palette = c("#cf3d3e", "#403990")) {
  data <- read.csv(input_csv)
  data <- data[, !colnames(data) %in% "sample_id"]
  data <- na.omit(data)

  cox_model <- coxph(Surv(time, status) ~ ., data = data)
  summary_cox <- summary(cox_model)

  data$risk_score <- predict(cox_model, type = "lp")
  data$risk_group <- ifelse(data$risk_score > median(data$risk_score), "high", "low")
  data$time <- data$time / 30  # convert days to months

  km_fit <- survfit(Surv(time, status) ~ risk_group, data = data)
  cox_result <- coxph(Surv(time, status) ~ risk_group, data = data)
  cox_sum <- summary(cox_result)

  HR_low <- round(cox_sum$coefficients[1, 2], 3)
  HR_high <- round(1 / HR_low, 3)
  CI_low <- round(1 / cox_sum$conf.int[1, 4], 3)
  CI_high <- round(1 / cox_sum$conf.int[1, 3], 3)
  p_val <- formatC(cox_sum$coefficients[1, 5], format = "e", digits = 2)

  svg(output_svg, width = 8, height = 6)
  km_plot <- ggsurvplot(
    km_fit, data = data, xlim = c(0, 120), break.x.by = 12,
    ylab = "Survival Probability", xlab = "Time (Months)",
    surv.median.line = "hv", legend.labs = c("High Risk Group", "Low Risk Group"),
    legend.title = "", surv.scale = "percent", palette = palette,
    title = title, risk.table = TRUE, risk.table.title = "Number at risk"
  )

  km_plot$plot <- km_plot$plot +
    annotate("text", x = 20, y = 0.1,
             label = paste0("HR = ", HR_high, " (95% CI: ", CI_low, "–", CI_high, ")\n", "p = ", p_val),
             size = 4, fontface = "plain")

  print(km_plot)
  dev.off()
}

#### 04. Generate KM Plots for All Pathways ####
plot_km_survival("COX_DATA.csv", "km_all.svg", "Metabolic Pathways")
plot_km_survival("lysine_cox.csv", "km_lysine.svg", "Lysine Degradation")
plot_km_survival("Alanine_cox.csv", "km_Alanine.svg", "Alanine, Aspartate, and Glutamate Metabolism")
plot_km_survival("AA_cox.csv", "km_aa.svg", "Biosynthesis of Amino Acids")
plot_km_survival("Nucleotide_cox.csv", "km_Nucleotide.svg", "Nucleotide Metabolism")
plot_km_survival("purine_cox.csv", "km_purine.svg", "Purine Metabolism")
plot_km_survival("galactose_cox.csv", "km_galactose.svg", "Galactose Metabolism")
plot_km_survival("inositol_cox.csv", "km_inositol.svg", "Inositol Phosphate Metabolism")
plot_km_survival("cofactors_cox.csv", "km_cofactors.svg", "Biosynthesis of Cofactors")
plot_km_survival("Oglycan_cox.csv", "km_O-glycan.svg", "Other Types of O-Glycan Biosynthesis")
plot_km_survival("Taurine_cox.csv", "km_taurine.svg", "Taurine and Hypotaurine Metabolism")
plot_km_survival("Oxidative_cox.csv", "km_oxidative.svg", "Oxidative Phosphorylation")
plot_km_survival("Phosphatidylinositol_cox.csv", "km_Phosphatidylinositol.svg", "Phosphatidylinositol Signaling System")
plot_km_survival("Carbon_cox.csv", "km_carbon.svg", "Carbon Metabolism")
