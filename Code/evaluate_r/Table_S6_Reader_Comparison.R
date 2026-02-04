# Script Name: Table_S6_Reader_Comparison.R
# Description: 
#   1. Compares multiple predictors against a Reference.
#   2. Performs McNemar's test (Sens/Spec), Bootstrap (PPV/NPV), and DeLong (AUC).
#   3. Generates summary rows for Junior, Senior, and All pathologists.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
required_packages <- c("pROC", "dplyr", "stringr", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(pROC)
library(dplyr)
library(stringr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure all columns are treated as numeric for analysis
raw_data <- raw_data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

# 2. Helper Functions
# ------------------------------------------------------------------------------
# Format proportion with 95% Confidence Interval
fmt_ci <- function(k, n) {
  if (n > 0) {
    res <- prop.test(k, n, conf.level = 0.95)
    sprintf("%.3f (%.3f-%.3f)", k/n, res$conf.int[1], res$conf.int[2])
  } else "NA"
}

# Format P-values for publication
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# 3. Core Analysis Engine
# ------------------------------------------------------------------------------
analyze_comparator <- function(target, ref_pred, comp_pred, col_name, is_ref = FALSE) {
  
  # --- Metric Calculation ---
  cm <- table(factor(target, levels = c(0, 1)), factor(comp_pred, levels = c(0, 1)))
  TP <- cm["1", "1"]; FN <- cm["1", "0"]
  TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  val_sens <- fmt_ci(TP, TP + FN)
  val_spec <- fmt_ci(TN, TN + FP)
  val_ppv  <- fmt_ci(TP, TP + FP)
  val_npv  <- fmt_ci(TN, TN + FN)
  
  roc_comp <- roc(target, comp_pred, quiet = TRUE)
  auc_val  <- as.numeric(auc(roc_comp))
  auc_ci   <- tryCatch(ci.auc(roc_comp, method = "delong"), error = function(e) c(NA, NA, NA))
  val_auc  <- if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA"
  
  # --- Statistical Comparison against Reference ---
  p_sens <- "Ref."; p_spec <- "Ref."; p_ppv <- "Ref."; p_npv <- "Ref."; p_auc <- "Ref."
  
  if (!is_ref) {
    # McNemar's Test for Sensitivity and Specificity
    idx_pos <- which(target == 1)
    if(length(idx_pos) > 0) {
      tbl_sens <- table(factor(ref_pred[idx_pos], levels=c(0,1)), factor(comp_pred[idx_pos], levels=c(0,1)))
      p_sens <- fmt_p(tryCatch(mcnemar.test(tbl_sens)$p.value, error=function(e) NA))
    }
    idx_neg <- which(target == 0)
    if(length(idx_neg) > 0) {
      tbl_spec <- table(factor(ref_pred[idx_neg], levels=c(0,1)), factor(comp_pred[idx_neg], levels=c(0,1)))
      p_spec <- fmt_p(tryCatch(mcnemar.test(tbl_spec)$p.value, error=function(e) NA))
    }
    
    # DeLong's Test for AUC Comparison
    roc_ref <- roc(target, ref_pred, quiet = TRUE)
    test_res <- tryCatch(roc.test(roc_ref, roc_comp, method = "delong"), error=function(e) NULL)
    if(!is.null(test_res)) p_auc <- fmt_p(test_res$p.value)
    
    # Bootstrap for PPV and NPV comparison (1000 iterations)
    calc_diff <- function(data, idx) {
      d <- data[idx, ]
      cm_r <- table(factor(d$t, 0:1), factor(d$r, 0:1))
      ppv_r <- cm_r[2,2]/sum(cm_r[,2]); npv_r <- cm_r[1,1]/sum(cm_r[,1])
      cm_c <- table(factor(d$t, 0:1), factor(d$c, 0:1))
      ppv_c <- cm_c[2,2]/sum(cm_c[,2]); npv_c <- cm_c[1,1]/sum(cm_c[,1])
      c(ppv_c - ppv_r, npv_c - npv_r)
    }
    
    set.seed(123)
    boot_dat <- data.frame(t=target, r=ref_pred, c=comp_pred)
    boot_diffs <- replicate(1000, calc_diff(boot_dat, sample(nrow(boot_dat), replace=TRUE)))
    
    get_p <- function(vec) {
      p <- 2 * min(mean(vec > 0, na.rm=TRUE), mean(vec < 0, na.rm=TRUE))
      fmt_p(min(p, 1.0))
    }
    p_ppv <- get_p(boot_diffs[1,]); p_npv <- get_p(boot_diffs[2,])
  }
  
  data.frame(Reader = col_name, 
             Sensitivity = val_sens, P_Sens = p_sens,
             Specificity = val_spec, P_Spec = p_spec,
             PPV = val_ppv, P_PPV = p_ppv,
             NPV = val_npv, P_NPV = p_npv,
             AUROC = val_auc, P_AUC = p_auc, stringsAsFactors = FALSE)
}

# 4. Main Table Generation
# ------------------------------------------------------------------------------
target_col <- raw_data[[1]]
ref_col    <- raw_data[[2]]
col_names  <- colnames(raw_data)

message("Processing individual columns...")
results_list <- lapply(2:ncol(raw_data), function(i) {
  analyze_comparator(target_col, ref_col, raw_data[[i]], col_names[i], is_ref = (i == 2))
})

final_table <- do.call(rbind, results_list)

# 5. Summary Rows (Averages)
# ------------------------------------------------------------------------------
# Extract numeric value from formatted string "0.xxx (0.xxx-0.xxx)"
get_val <- function(str_val) {
  if (str_val == "NA" || is.na(str_val)) return(NA)
  as.numeric(str_split(str_val, " ")[[1]][1])
}

# Calculate average row based on specified indices
calc_avg_row <- function(df, row_indices, label_name) {
  metric_cols <- c("Sensitivity", "Specificity", "PPV", "NPV", "AUROC")
  new_row <- list(Reader = label_name)
  
  for (col in colnames(df)) {
    if (col == "Reader") next
    if (col %in% metric_cols) {
      vals <- sapply(df[row_indices, col], get_val)
      new_row[[col]] <- sprintf("%.3f", mean(vals, na.rm = TRUE))
    } else {
      new_row[[col]] <- "-" # No statistical test for average rows
    }
  }
  return(as.data.frame(new_row, stringsAsFactors = FALSE))
}

# Define grouping indices (Column 3-6: Junior, Column 7-9: Senior)
idx_junior <- (3:6) - 1
idx_senior <- (7:9) - 1
idx_all    <- (3:9) - 1

avg_junior <- calc_avg_row(final_table, idx_junior, "Average (Junior pathologists)")
avg_senior <- calc_avg_row(final_table, idx_senior, "Average (Senior pathologists)")
avg_all    <- calc_avg_row(final_table, idx_all,    "Average (All pathologists)")

final_table_with_avg <- rbind(final_table, avg_junior, avg_senior, avg_all)

# 6. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_table_with_avg, file.path(results_dir, "Table_S6_Comparison_Results.csv"), row.names = FALSE)

print(final_table_with_avg)
message("Analysis complete. Results exported to /results/.")


# ==============================================================================
