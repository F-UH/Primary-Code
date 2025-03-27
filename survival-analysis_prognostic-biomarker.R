### 0. Load Packages ====
packages <- c("TCGAbiolinks", "survminer", "survival", "SummarizedExperiment", 
              "tidyverse", "DESeq2", "stringr", "dbplyr", "ggsurvfit", "broom", "forestplot")
invisible(lapply(packages, library, character.only = TRUE))

### 1. Clinical Data Processing ====
clinical_luad <- GDCquery_clinic("TCGA-LUAD")
clinical_lusc <- GDCquery_clinic("TCGA-LUSC")

# Create unified survival variables
for (df in list(clinical_luad, clinical_lusc)) {
  df$deceased <- ifelse(df$vital_status == "Alive", FALSE, TRUE)
  df$overall_survival <- ifelse(df$deceased, df$days_to_death, df$days_to_last_follow_up)
}

# Extract useful columns and rename
colnames_target <- c("sample_id", "ajcc_pathologic_stage", "subtype", "race", "sex", "ethnicity",
                     "vital_status", "deceased", "survival_time", "age", "smoking_status")

clinical_luad_extract <- clinical_luad %>%
  transmute(sample_id = submitter_id, ajcc_pathologic_stage, subtype = primary_diagnosis,
            race, sex = gender, ethnicity, vital_status, deceased, 
            survival_time = overall_survival, age = age_at_index, 
            smoking_status = pack_years_smoked, primary_diagnosis = "LUAD")

clinical_lusc_extract <- clinical_lusc %>%
  transmute(sample_id = submitter_id, ajcc_pathologic_stage, subtype = primary_diagnosis,
            race, sex = gender, ethnicity, vital_status, deceased, 
            survival_time = overall_survival, age = age_at_index, 
            smoking_status = pack_years_smoked, primary_diagnosis = "LUSC")

nsclc_clinical <- bind_rows(clinical_luad_extract, clinical_lusc_extract)
nsclc_clinical$sample_id <- gsub("-", "_", nsclc_clinical$sample_id)

### 2. TPM tRNA Expression Matrix ====
tpm <- read.csv("TCGA_tRNA_tpms.csv", header = TRUE, row.names = 1)
tpm <- tpm[rowSums(tpm == 0) / ncol(tpm) < 0.5, ]

# Select 6 tRNAs
target_genes <- c("tRNA-Val-CAC-2-1", "tRNA-Leu-AAG-2-3", "tRNA-Lys-CTT-3-1", 
                  "tRNA-Val-CAC-1-5", "tRNA-Ala-TGC-3-2", "tRNA-Asp-GTC-1-1")
tpm_filtered <- tpm[rownames(tpm) %in% target_genes, ]

colnames(tpm_filtered) <- str_extract(colnames(tpm_filtered), "^[^_]+_[^_]+_[^_]+")
matched_clinical <- nsclc_clinical[nsclc_clinical$sample_id %in% colnames(tpm_filtered), ]
tpm_long <- as.data.frame(t(tpm_filtered))
tpm_long$sample_id <- rownames(tpm_long)

meta_data <- merge(tpm_long, matched_clinical, by = "sample_id")

### 3. KM Plot Function ====
plot_single_km <- function(file, cutoff, gene_name, color = c("#403990", "#cf3d3e")) {
  df <- read.csv(file)
  df$compare <- ifelse(df$compare >= cutoff, 1, 0)
  df$time <- df$time / 30

  svg(paste0(gene_name, "_plot.svg"), width = 8, height = 6)
  km <- survfit(Surv(time, event) ~ compare, data = df)
  ggsurvplot(km, data = df, xlim = c(0, 72), break.x.by = 12,
             ylab = "Survival Probability", xlab = "Time (Months)",
             pval = TRUE, surv.median.line = "hv",
             legend.labs = c(paste0("Low ", gene_name), paste0("High ", gene_name)),
             legend.title = "", palette = color, title = "",
             risk.table = TRUE, risk.table.title = "Number at risk")
  dev.off()
}

# Draw KM plots
plot_single_km("tRNA-Lys-CTT-3-1.csv", 561, "tRNA-Lys-CTT-3-1")
plot_single_km("tRNA-Val-CAC-2-1.csv", 370, "tRNA-Val-CAC-2-1")
plot_single_km("tRNA-Leu-AAG-2-3.csv", 4492, "tRNA-Leu-AAG-2-3")

### 4. Combine Signature Score ====
meta_data$`tRNA-Leu-AAG-2-3` <- ifelse(meta_data$`tRNA-Leu-AAG-2-3` >= 4492, 1, 0)
meta_data$`tRNA-Val-CAC-2-1` <- ifelse(meta_data$`tRNA-Val-CAC-2-1` >= 370, 1, 0)
meta_data$`tRNA-Lys-CTT-3-1` <- ifelse(meta_data$`tRNA-Lys-CTT-3-1` >= 561, 1, 0)
meta_data$deceased <- as.numeric(meta_data$deceased)
meta_data$survival_time <- meta_data$survival_time / 30

meta_data$signature_score <- rowSums(meta_data[, c("tRNA-Leu-AAG-2-3", "tRNA-Val-CAC-2-1", "tRNA-Lys-CTT-3-1")])

# Optimal cutpoint
opt_cut <- surv_cutpoint(meta_data, time = "survival_time", event = "deceased", variables = "signature_score")
print(opt_cut)

meta_data$risk_group <- ifelse(meta_data$signature_score > 1, "High", "Low")
svg("combine_tRNA.svg", width = 8, height = 6)
fit <- survfit(Surv(survival_time, deceased) ~ risk_group, data = meta_data)
ggsurvplot(fit, xlim = c(0, 72), break.x.by = 12, pval = TRUE, surv.median.line = "hv",
           ylab = "Survival Probability", xlab = "Time (Months)",
           legend.labs = c("High Risk Score", "Low Risk Score"),
           legend.title = "", palette = c("#fcbb44", "#839dd1"),
           title = "", risk.table = TRUE, risk.table.title = "Number at risk")
dev.off()

### 5. Cox Regression Model and Forest Plot ====
meta_data$smoking_status <- ifelse(is.na(meta_data$smoking_status), 0, 1)
meta_data$risk_group <- relevel(factor(meta_data$risk_group), ref = "Low")
meta_data$sex <- relevel(factor(meta_data$sex), ref = "female")
meta_data$ajcc_pathologic_stage <- relevel(factor(meta_data$ajcc_pathologic_stage), ref = "Stage I")

cox_model <- coxph(Surv(survival_time, deceased) ~ risk_group + subtype + sex + age + 
                     smoking_status + ajcc_pathologic_stage, data = meta_data)

tidied <- broom::tidy(cox_model)
tidied$HR <- exp(tidied$estimate)
tidied$lower_CI <- exp(tidied$estimate - 1.96 * tidied$std.error)
tidied$upper_CI <- exp(tidied$estimate + 1.96 * tidied$std.error)

# Manually construct display table
label_table <- cbind(
  c("Subgroup", tidied$term),
  c("N", rep("", nrow(tidied))),  # placeholder
  c("Adjusted HR (95% CI)", sprintf("%.2f (%.2f–%.2f)", tidied$HR, tidied$lower_CI, tidied$upper_CI)),
  c("p-value", signif(tidied$p.value, 3)),
  c("Significance", ifelse(tidied$p.value < 0.05, "*", ""))
)

svg("forest_plot.svg", width = 10, height = 6)
forestplot(labeltext = label_table,
           mean = c(NA, tidied$HR),
           lower = c(NA, tidied$lower_CI),
           upper = c(NA, tidied$upper_CI),
           zero = 1, boxsize = 0.3, lineheight = unit(6, "mm"),
           graph.pos = 4, colgap = unit(5, "mm"),
           xticks = c(0.5, 1, 1.5),
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.8), cex = 0.9),
           col = fpColors(box = "#110075", summary = "black", lines = "black", zero = "grey60"),
           hrzl_lines = list("1" = gpar(lty = 1), "2" = gpar(lty = 1), nrow(label_table) = gpar(lty = 1)))
dev.off()
