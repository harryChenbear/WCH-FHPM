# Script Name: Table_1_Performance_Analysis.R
# Description: Calculate comprehensive performance metrics (Sens, Spec, PPV, 
#              NPV, Brier Score, AUROC) including an "Overall" dataset.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)   # For ROC/AUC analysis
library(dplyr)  # For data manipulation
library(here)   # For relative path management

# Load data using relative path
data_path <- here("data", "Table 1.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preparation
# ------------------------------------------------------------------------------
# Create "Overall Dataset" by duplicating all entries
overall_data <- raw_data
overall_data$Dataset <- "Overall Dataset"

# Combine individual datasets with the overall summary
combined_data <- rbind(raw_data, overall_data)

# Set factor levels to control the order (Individual datasets followed by Overall)
orig_names <- unique(raw_data$Dataset)
combined_data$Dataset <- factor(combined_data$Dataset, levels = c(orig_names, "Overall Dataset"))

# 3. Metric Calculation Function
# ------------------------------------------------------------------------------
calculate_metrics_final <- function(group_data) {
  
  # Ensure labels and predictions are numeric
  actual <- as.numeric(as.character(group_data$true_labels))
  predicted <- as.numeric(as.character(group_data$pred_labels)) 
  
  # --- Part 1: Sensitivity, Specificity, PPV, NPV ---
  # Construct confusion matrix
  cm <- table(factor(actual, levels = c(0, 1)), 
              factor(predicted, levels = c(0, 1)))
  
  TP <- cm["1", "1"]; FN <- cm["1", "0"]
  TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  # Helper for calculating proportions with 95% CI
  get_prop_ci <- function(k, n) {
    if (n > 0) {
      res <- prop.test(k, n, conf.level = 0.95)
      return(c(val = k/n, lower = res$conf.int[1], upper = res$conf.int[2]))
    } else {
      return(c(val = 0, lower = 0, upper = 0))
    }
  }
  
  sens <- get_prop_ci(TP, TP + FN)
  spec <- get_prop_ci(TN, TN + FP)
  ppv  <- get_prop_ci(TP, TP + FP)
  npv  <- get_prop_ci(TN, TN + FN)
  
  # --- Part 2: AUROC (DeLong Method) ---
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- ci.auc(roc_obj, conf.level = 0.95, method = "delong")
  
  # --- Part 3: Brier Score (Bootstrap CI) ---
  brier_val <- mean((predicted - actual)^2)
  
  set.seed(123) # Ensure reproducibility for bootstrap
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), length(actual), replace = TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs = c(0.025, 0.975))
  
  # --- Part 4: Formatting Result ---
  fmt <- function(val, low, high) {
    sprintf("%.3f (%.3f-%.3f)", val, low, high)
  }
  
  result <- data.frame(
    Dataset = unique(group_data$Dataset),
    Sensitivity = fmt(sens['val'], sens['lower'], sens['upper']),
    Specificity = fmt(spec['val'], spec['lower'], spec['upper']),
    PPV         = fmt(ppv['val'],  ppv['lower'],  ppv['upper']),
    NPV         = fmt(npv['val'],  npv['lower'],  npv['upper']),
    Brier_Score = fmt(brier_val, brier_ci[1], brier_ci[2]),
    AUROC       = fmt(auc_val, auc_ci[1], auc_ci[3]),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# 4. Execution and Export
# ------------------------------------------------------------------------------
final_results <- combined_data %>%
  group_by(Dataset) %>%
  do(calculate_metrics_final(.)) %>%
  ungroup()

# Export results to the 'results' directory
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_results, file.path(results_dir, "Table_1_Performance_Metrics.csv"), row.names = FALSE)

# Display results in console
print(final_results)
message("Calculation complete. Results exported to /results/ directory.")


# ==============================================================================
