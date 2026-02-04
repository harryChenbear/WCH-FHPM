# Script Name: Table_S7_Paired_Analysis.R
# Description: Paired performance analysis by pathologist group.
#              Compares "Modified strategy" vs "Without AI".
#              Metrics: Sens, Spec, PPV, NPV, AUROC, Brier Score.
#              Tests: McNemar, DeLong, and Bootstrap for paired differences.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(stringr)
library(here)

# Load data using relative path
# check.names = FALSE ensures special characters in headers are preserved
data_path <- here("data", "Table S7.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path, check.names = FALSE)

# Ensure numeric format across all diagnostic columns
raw_data <- raw_data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

# 2. Helper Functions for Formatting and Stats
# ------------------------------------------------------------------------------
fmt_ci <- function(k, n) {
  if (n > 0) {
    res <- prop.test(k, n, conf.level = 0.95)
    sprintf("%.3f (%.3f-%.3f)", k/n, res$conf.int[1], res$conf.int[2])
  } else "NA"
}

fmt_val_ci <- function(val, lower, upper) {
  sprintf("%.3f (%.3f-%.3f)", val, lower, upper)
}

fmt_p <- function(p) {
  if (is.na(p) || is.nan(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# 3. Paired Group Analysis Function
# ------------------------------------------------------------------------------
analyze_paired_group <- function(target, ref_pred, comp_pred, group_name) {
  
  # --- A. Metrics for "Without AI" (Reference) ---
  cm_r <- table(factor(target, levels=0:1), factor(ref_pred, levels=0:1))
  sens_r <- fmt_ci(cm_r[2,2], sum(cm_r[2,])); spec_r <- fmt_ci(cm_r[1,1], sum(cm_r[1,]))
  ppv_r  <- fmt_ci(cm_r[2,2], sum(cm_r[,2])); npv_r  <- fmt_ci(cm_r[1,1], sum(cm_r[,1]))
  
  roc_r <- roc(target, ref_pred, quiet=TRUE)
  auc_r_val <- as.numeric(auc(roc_r))
  auc_r_ci  <- tryCatch(ci.auc(roc_r, method="delong"), error=function(e) c(NA,NA,NA))
  auc_r <- fmt_val_ci(auc_r_val, auc_r_ci[1], auc_r_ci[3])
  
  brier_r_val <- mean((ref_pred - target)^2)
  set.seed(123)
  b_boot_r <- replicate(1000, { idx<-sample(length(target), replace=T); mean((ref_pred[idx]-target[idx])^2) })
  brier_r_ci <- quantile(b_boot_r, c(0.025, 0.975))
  brier_r <- fmt_val_ci(brier_r_val, brier_r_ci[1], brier_r_ci[2])
  
  # --- B. Metrics for "Modified strategy" (Comparator) ---
  cm_c <- table(factor(target, levels=0:1), factor(comp_pred, levels=0:1))
  sens_c <- fmt_ci(cm_c[2,2], sum(cm_c[2,])); spec_c <- fmt_ci(cm_c[1,1], sum(cm_c[1,]))
  ppv_c  <- fmt_ci(cm_c[2,2], sum(cm_c[,2])); npv_c  <- fmt_ci(cm_c[1,1], sum(cm_c[,1]))
  
  roc_c <- roc(target, comp_pred, quiet=TRUE)
  auc_c_val <- as.numeric(auc(roc_c))
  auc_c_ci  <- tryCatch(ci.auc(roc_c, method="delong"), error=function(e) c(NA,NA,NA))
  auc_c <- fmt_val_ci(auc_c_val, auc_c_ci[1], auc_c_ci[3])
  
  brier_c_val <- mean((comp_pred - target)^2)
  set.seed(123)
  b_boot_c <- replicate(1000, { idx<-sample(length(target), replace=T); mean((comp_pred[idx]-target[idx])^2) })
  brier_c_ci <- quantile(b_boot_c, c(0.025, 0.975))
  brier_c <- fmt_val_ci(brier_c_val, brier_c_ci[1], brier_c_ci[2])
  
  # --- C. Statistical Comparisons ---
  idx_pos <- which(target==1); idx_neg <- which(target==0)
  
  p_sens <- tryCatch(mcnemar.test(table(factor(ref_pred[idx_pos],0:1), factor(comp_pred[idx_pos],0:1)))$p.value, error=function(e) NA)
  p_spec <- tryCatch(mcnemar.test(table(factor(ref_pred[idx_neg],0:1), factor(comp_pred[idx_neg],0:1)))$p.value, error=function(e) NA)
  p_auc <- tryCatch(roc.test(roc_r, roc_c, method="delong")$p.value, error=function(e) NA)
  
  # Bootstrap for PPV, NPV, and Brier Score differences
  calc_diffs <- function(d, idx) {
    curr <- d[idx, ]
    tb_r <- table(factor(curr$t,0:1), factor(curr$r,0:1))
    ppv_r <- tb_r[2,2]/sum(tb_r[,2]); npv_r <- tb_r[1,1]/sum(tb_r[,1]); br_r <- mean((curr$r - curr$t)^2)
    tb_c <- table(factor(curr$t,0:1), factor(curr$c,0:1))
    ppv_c <- tb_c[2,2]/sum(tb_c[,2]); npv_c <- tb_c[1,1]/sum(tb_c[,1]); br_c <- mean((curr$c - curr$t)^2)
    c(ppv_c - ppv_r, npv_c - npv_r, br_c - br_r)
  }
  
  boot_df <- data.frame(t=target, r=ref_pred, c=comp_pred)
  set.seed(123)
  boot_res <- replicate(1000, calc_diffs(boot_df, sample(nrow(boot_df), replace=TRUE)))
  
  get_p <- function(vec) {
    p <- 2 * min(mean(vec > 0, na.rm=TRUE), mean(vec < 0, na.rm=TRUE))
    fmt_p(min(p, 1.0))
  }
  p_ppv <- get_p(boot_res[1,]); p_npv <- get_p(boot_res[2,]); p_brier <- get_p(boot_res[3,])
  
  # --- D. Formatting Output Table ---
  row_comp <- data.frame(
    Group = group_name, Strategy = "Modified strategy",
    Sensitivity = sens_c, P_Sens = fmt_p(p_sens),
    Specificity = spec_c, P_Spec = fmt_p(p_spec),
    PPV = ppv_c, P_PPV = p_ppv,
    NPV = npv_c, P_NPV = p_npv,
    AUROC = auc_c, P_AUC = fmt_p(p_auc),
    Brier_Score = brier_c, P_Brier = p_brier,
    stringsAsFactors = FALSE
  )
  
  row_ref <- data.frame(
    Group = group_name, Strategy = "Without AI",
    Sensitivity = sens_r, P_Sens = "Ref.",
    Specificity = spec_r, P_Spec = "Ref.",
    PPV = ppv_r, P_PPV = "Ref.",
    NPV = npv_r, P_NPV = "Ref.",
    AUROC = auc_r, P_AUC = "Ref.",
    Brier_Score = brier_r, P_Brier = "Ref.",
    stringsAsFactors = FALSE
  )
  
  return(rbind(row_comp, row_ref))
}

# 4. Processing Loop
# ------------------------------------------------------------------------------
target_col <- raw_data[[1]]
all_colnames <- colnames(raw_data)[-1]
groups <- unique(sapply(str_split(all_colnames, "_"), `[`, 1))

results_list <- list()

message("Executing paired analysis...")

for (grp in groups) {
  cols_in_grp <- all_colnames[startsWith(all_colnames, paste0(grp, "_"))]
  col_ref <- cols_in_grp[grep("Without AI", cols_in_grp)]
  col_comp <- cols_in_grp[grep("Modified strategy", cols_in_grp)]
  
  if (length(col_ref) >= 1 && length(col_comp) >= 1) {
    results_list[[grp]] <- analyze_paired_group(target_col, raw_data[[col_ref[1]]], raw_data[[col_comp[1]]], grp)
  }
}

final_output <- do.call(rbind, results_list)

# 5. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_output, file.path(results_dir, "Table_S7_Paired_Comparison.csv"), row.names = FALSE)

print(final_output)
message("Analysis complete. Results exported to /results/.")


# ==============================================================================
