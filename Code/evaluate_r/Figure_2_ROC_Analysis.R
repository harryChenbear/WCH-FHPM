# ==============================================================================
# Script Name: Figure_2_ROC_Analysis.R
# Description: Generate ROC curves for multiple datasets.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
# Automatically install and load required packages
required_packages <- c("ggplot2", "pROC", "dplyr", "gridExtra", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here) 

# 2. Data Import
# ------------------------------------------------------------------------------
# Load data from the 'data' directory
data_path <- here("data", "Figure 2.csv")

if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 2.csv' is in the /data/ directory.")
}

raw_data <- read.csv(data_path)

# 3. Data Preprocessing
# ------------------------------------------------------------------------------
datasets <- unique(raw_data$Dataset)
raw_data$Dataset <- factor(raw_data$Dataset, levels = datasets)

# Define color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 4. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list() 

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (AUC & CI) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  # Construct annotation string
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Coordinate Extraction ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  # --- Step C: Create Plot Object ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = custom_colors[ds], linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             size = 4, hjust = 0, lineheight = 1.2) +
    
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = ds) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[ds]] <- p
}

# 5. Export Results
# ------------------------------------------------------------------------------
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(output_dir, "Figure2_ROC_Analysis.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Output saved to: ", output_dir)




# ==============================================================================
