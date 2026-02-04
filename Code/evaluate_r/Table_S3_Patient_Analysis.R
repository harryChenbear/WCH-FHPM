# Script Name: Table_S3_Patient_Analysis.R
# Description: Calculate aggregated patient-level performance metrics.
#              Rule: Max pred_label per patient is used for diagnosis.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S3.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure labels are numeric
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# 2. Patient-Level Aggregation
# ------------------------------------------------------------------------------
# Robust identification of Patient ID column (matching "Pts")
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (is.na(pts_col)) stop("Patient ID column (Pts) not found.")
colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Aggregation: Assign patient-level diagnosis based on the maximum label per patient
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>% 
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Add Overall Summary
# ------------------------------------------------------------------------------
overall_data <- patient_data
overall_data$Dataset <- "Overall Dataset"

combined_data <- rbind(patient_data, overall_data)

# Control dataset order: Individual sets first, followed by Overall summary
orig_names <- unique(patient_data$Dataset)
combined_data$Dataset <- factor(combined_data$Dataset, levels = c(orig_names, "Overall Dataset"))

# 4. Metrics Calculation Function
# ------------------------------------------------------------------------------
calc_metrics <- function(df) {
  
  actual <- df$true_labels
  predicted <- df$pred_labels
  
  # --- Contingency Table Metrics (Sens, Spec, PPV, NPV) ---
  cm <- table(factor(actual, levels=c(0,1)), factor(predicted, levels=c(0,1)))
  TP <- cm["1","1"]; FN <- cm["1","0"]; TN <- cm["0","0"]; FP <- cm["0","1"]
  
  # Proportion test for confidence intervals
  get_ci <- function(k, n) {
    if(n > 0) {
      res <- prop.test(k, n, conf.level=0.95)
      c(res$estimate, res$conf.int)
    } else c(0, 0, 0)
  }
  
  sens <- get_ci(TP, TP+FN)
  spec <- get_ci(TN, TN+FP)
  ppv  <- get_ci(TP, TP+FP)
  npv  <- get_ci(TN, TN+FN)
  
  # --- AUROC (DeLong Method) ---
  roc_obj <- roc(actual, predicted, quiet=TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  # Using tryCatch to handle edge cases where DeLong might fail
  auc_ci <- tryCatch(ci.auc(roc_obj, method="delong"), error=function(e) c(NA,NA,NA))
  
  # --- Brier Score (Bootstrap) ---
  brier_val <- mean((predicted - actual)^2)
  set.seed(123) # Reproducibility
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), replace=TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs=c(0.025, 0.975))
  
  # --- Formatting Results ---
  fmt <- function(v) sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
  
  data.frame(
    Dataset = unique(df$Dataset),
    Sensitivity = fmt(sens),
    Specificity = fmt(spec),
    PPV = fmt(ppv),
    NPV = fmt(npv),
    Brier_Score = sprintf("%.3f (%.3f-%.3f)", brier_val, brier_ci[1], brier_ci[2]),
    AUROC = if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA",
    stringsAsFactors = FALSE
  )
}

# 5. Execution and Export
# ------------------------------------------------------------------------------
final_results <- combined_data %>%
  group_by(Dataset) %>%
  do(calc_metrics(.)) %>%
  ungroup()

# Create results directory if needed
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_results, file.path(results_dir, "Table_S3_Patient_Metrics.csv"), row.names = FALSE)

# Display Results
print(final_results)
message("Patient-level metrics calculation complete. Results saved to /results/.")


# ==============================================================================
