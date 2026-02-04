# Script Name: Figure_S6_Patient_ROC.R
# Description: Generate Patient-Level ROC curves for Figure S6.
#              - Aggregation: Maximum probability and label per Patient ID.
#              - Curve Shape: Derived from 'pred_probs'.
#              - AUC Statistics: Derived from 'pred_labels'.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
# Robustly identify Patient ID column
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (!is.na(pts_col)) colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Ensure numeric types
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_probs  <- as.numeric(as.character(raw_data$pred_probs))

# Identify prediction labels column (defaults to 5th column if name not found)
if(!"pred_labels" %in% colnames(raw_data)){
  colnames(raw_data)[5] <- "pred_labels"
}
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# Patient-Level Aggregation: Taking the maximum values
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>%
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_probs  = max(pred_probs, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

datasets <- unique(patient_data$Dataset)
patient_data$Dataset <- factor(patient_data$Dataset, levels = datasets)

# Define color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 3. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(patient_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (Based on binary labels) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Coordinate Extraction (Based on probabilities for smooth curve) ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  current_color <- custom_colors[as.character(ds)]
  if(is.null(current_color)) current_color <- "black"
  
  # --- Step C: Create Individual Plot ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = current_color, linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 4, hjust = 0, lineheight = 1.2) +
    
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
  
  plot_list[[as.character(ds)]] <- p
}

# 4. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure_S6_Patient_ROC.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Results exported to /results/.")


# ==============================================================================
