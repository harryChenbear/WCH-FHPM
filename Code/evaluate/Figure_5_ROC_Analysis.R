# Script Name: Figure_5_ROC_Analysis.R
# Description: Generate ROC curves for Figure 5.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here)

# Load data from the 'data' directory
data_path <- here("data", "Figure 5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
datasets <- unique(raw_data$Dataset)
raw_data$Dataset <- factor(raw_data$Dataset, levels = datasets)

# Define consistent color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 3. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list() 

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (AUC & CI) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  # Prepare annotation text
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Generate ROC Coordinates ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  # --- Step C: Create Individual Plot ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    # Draw curve with standardized linewidth
    geom_path(color = custom_colors[ds], linewidth = 0.8) +
    
    # Reference diagonal line
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 4, hjust = 0, lineheight = 1.2) +
    
    # Configure axes with reversed specificity and padding
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

# 4. Final Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure5_ROC_Analysis.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Combined plot saved to: ", results_dir)


# ==============================================================================
