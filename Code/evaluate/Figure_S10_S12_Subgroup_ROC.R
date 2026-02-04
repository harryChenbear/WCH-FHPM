# Script Name: Figure_S10_S12_Subgroup_ROC.R
# Description: Generate ROC curves for subgroups in Figure S10 and S12.
#              - Curve Shape: Based on 'pred_probs' (Continuous).
#              - Statistics: AUROC & 95% CI based on 'pred_labels' (Binary).
#              - Layout: Automated grid arrangement based on group count.
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
data_path <- here("data", "Figure S10 S12.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure numeric format for columns
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))
raw_data$pred_probs  <- as.numeric(as.character(raw_data$pred_probs))

groups <- unique(raw_data$GROUP)

# 2. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list()
line_color <- "#1C93B5" 

for (grp in groups) {
  sub_data <- subset(raw_data, GROUP == grp)
  
  # --- Part A: Statistical Computation (Binary Labels) ---
  # AUC calculated using binary prediction labels as requested
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  
  ci_result <- tryCatch({
    ci.auc(roc_stat_obj, method = "delong")
  }, error = function(e) {
    return(c(NA, NA, NA))
  })
  
  # Format annotation text
  if (all(is.na(ci_result))) {
    label_str <- "AUC: NA"
  } else {
    label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                        "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  }
  
  # Define text position df
  anno_df <- data.frame(x = 0.5, y = 0.15, label = label_str)
  
  # --- Part B: Coordinate Extraction (Continuous Probabilities) ---
  # ROC curve shape based on predicted probabilities
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Group = grp
  )
  
  # --- Part C: Individual Plot Creation ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = line_color, linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 3.5, hjust = 0, lineheight = 1.2) +
    
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = grp) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[grp]] <- p
}

# 3. Grid Assembly and Export
# ------------------------------------------------------------------------------
# Calculate grid dimensions automatically
n_plots <- length(plot_list)
ncol_grid <- ceiling(sqrt(n_plots))
nrow_grid <- ceiling(n_plots / ncol_grid)

combined_plot <- grid.arrange(grobs = plot_list, ncol = ncol_grid)

# Define output path
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(filename = file.path(results_dir, "Figure_S10_S12_ROC.pdf"), 
       plot = combined_plot, 
       width = ncol_grid * 3.5, 
       height = nrow_grid * 3.5, 
       limitsize = FALSE, device = "pdf")

message("Subgroup ROC analysis complete. Plot saved to /results/.")


# ==============================================================================
