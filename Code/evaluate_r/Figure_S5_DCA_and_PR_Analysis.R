# Script Name: Figure_S5_DCA_and_PR_Analysis.R
# Description: Generates Decision Curves and Precision-Recall Curves 
#              across multiple cohorts to evaluate clinical utility.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggthemes)
library(rmda)
library(pROC)
library(PRROC)
library(here)

# Define dataset and output paths
data_path <- here("data", "Figure S5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# 2. Configuration
# ------------------------------------------------------------------------------
group_col <- "Dataset" 
my_colors <- c("#7AB5D1", "#678655", "#F0A72C")

# Load Data
demo_data <- read.table(data_path, sep=",", header = TRUE)
if (!group_col %in% names(demo_data)) {
  stop(paste0("Error: Column '", group_col, "' not found."))
}

cohorts <- unique(demo_data[[group_col]])

# 3. PDF Export Initialization
# ------------------------------------------------------------------------------
output_filename <- file.path(results_dir, "Figure_S5_DCA_PR_Analysis.pdf")
message("Generating combined analysis PDF...")

pdf(output_filename, width = 10, height = 15) 

# Layout: 3 Rows, 2 Columns with outer margins
par(mfrow = c(3, 2), oma = c(2, 4, 1, 1))

# 4. Processing Loop
# ------------------------------------------------------------------------------
for (i in seq_along(cohorts)) {
  
  cohort_name <- cohorts[i]
  current_color <- my_colors[(i - 1) %% length(my_colors) + 1]
  sub_data <- demo_data[demo_data[[group_col]] == cohort_name, ]
  
  # --- Part A: Decision Curve Analysis (DCA) ---
  # Adjust margins for axis labels
  par(mar = c(7, 5, 3, 2)) 
  
  dca_model <- decision_curve(formula = true_labels ~ pred_probs, 
                              data = sub_data, family = binomial, 
                              thresholds = seq(0, 1, by = 0.01), policy = "opt-in")
  dca_model <- na.omit(dca_model)
  
  plot_decision_curve(dca_model,
                      curve.names = "Prediction Model",
                      col = current_color,
                      lwd = 3, lty = 1,
                      xlab = "Threshold Probability",
                      ylab = "Net Benefit",
                      cex.lab = 1.6, cex.axis = 1.4, cex.main = 1.8,
                      legend.position = "none",
                      cost.benefit.axis = TRUE)
  
  title("Decision Curve", col.main = current_color, font.main = 4, line = 1, cex.main = 1.8)
  
  # Add Dataset Label to the left margin
  mtext(cohort_name, side = 2, line = 5, cex = 1.5, font = 2, col = "black", las = 0)
  
  # --- Part B: Precision-Recall (PR) Curve ---
  par(mar = c(7, 5, 3, 2)) 
  
  true_labels <- sub_data$true_labels
  pred_probs <- sub_data$pred_probs
  
  pr <- pr.curve(scores.class0 = pred_probs[true_labels == 1],
                 scores.class1 = pred_probs[true_labels == 0], curve = TRUE)
  
  # Calculate ROC statistics for annotation
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE)
  auc_val <- auc(roc_obj)
  ci_val <- ci.auc(roc_obj)
  
  plot(pr,
       xlab = "Recall",
       ylab = "Precision",
       col = current_color,
       lwd = 3,
       cex.lab = 1.6, cex.axis = 1.4, cex.main = 1.8)
  
  # Annotate with AUC and Confidence Intervals
  text(x = 0, y = 0.05, 
       labels = sprintf("AUC: %.3f (95%% CI: %.3f-%.3f)", auc_val, ci_val[1], ci_val[3]),
       col = "black", cex = 1.3, adj = 0)
  
  box(lwd = 2)
  title("Recall Curve", col.main = current_color, font.main = 4, line = 1, cex.main = 1.8)
}

dev.off()
message("Analysis complete. PDF saved to: ", output_filename)


# ==============================================================================
