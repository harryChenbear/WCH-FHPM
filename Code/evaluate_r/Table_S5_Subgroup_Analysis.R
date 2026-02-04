# Script Name: Table_S5_Subgroup_Analysis.R
# Description: Subgroup performance metrics for Table S5.
#              Calculates Sens, Spec, PPV, NPV, Brier Score, and AUROC
#              using binary 'pred_labels'.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure numeric format for labels
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# 2. Metric Calculation Function
# ------------------------------------------------------------------------------
calculate_metrics <- function(df) {
  
  actual    <- df$true_labels
  predicted <- df$pred_labels 
  
  # --- Contingency Table Metrics (Sens, Spec, PPV, NPV) ---
  cm <- table(factor(actual, levels = c(0, 1)), factor(predicted, levels = c(0, 1)))
  TP <- cm["1", "1"]; FN <- cm["1", "0"]; TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  # Proportion test for 95% Confidence Intervals
  get_ci <- function(k, n) {
    if (n > 0) {
      res <- prop.test(k, n, conf.level = 0.95)
      c(res$estimate, res$conf.int)
    } else c(NA, NA, NA)
  }
  
  sens <- get_ci(TP, TP + FN)
  spec <- get_ci(TN, TN + FP)
  ppv  <- get_ci(TP, TP + FP)
  npv  <- get_ci(TN, TN + FN)
  
  # --- AUROC (DeLong Method using Binary Labels) ---
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) c(NA, NA, NA))
  
  # --- Brier Score (Bootstrap) ---
  brier_val <- mean((predicted - actual)^2)
  
  set.seed(123) # For reproducibility
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), replace = TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs = c(0.025, 0.975))
  
  # --- Formatting Results ---
  fmt <- function(v) {
    if(all(is.na(v))) return("NA")
    sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
  }
  
  data.frame(
    Group       = unique(df$GROUP),
    Sensitivity = fmt(sens),
    Specificity = fmt(spec),
    PPV         = fmt(ppv),
    NPV         = fmt(npv),
    Brier_Score = sprintf("%.3f (%.3f-%.3f)", brier_val, brier_ci[1], brier_ci[2]),
    AUROC       = if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA",
    stringsAsFactors = FALSE
  )
}

# 3. Execution and Export
# ------------------------------------------------------------------------------
final_results <- raw_data %>%
  group_by(GROUP) %>%
  do(calculate_metrics(.)) %>%
  ungroup()

# Create results directory if needed
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Save to CSV
write.csv(final_results, file.path(results_dir, "Table_S5_Metrics_Binary.csv"), row.names = FALSE)

# Display Results
print(final_results)
message("Subgroup analysis complete. Results saved to /results/.")


# ==============================================================================
