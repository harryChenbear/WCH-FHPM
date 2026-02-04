# Script Name: Figure_5_Confusion_Matrices.R
# Description: Generate confusion matrices for datasets in Figure 5.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "caret", "gridExtra", "scales", "dplyr", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(caret)
library(gridExtra)
library(scales)
library(dplyr)
library(here)

# 2. Data Import
# ------------------------------------------------------------------------------
# Load data from the 'data' directory
data_path <- here("data", "Figure 5.csv")

if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 5.csv' is in the /data/ directory.")
}

raw_data <- read.csv(data_path)

# 3. Data Preprocessing
# ------------------------------------------------------------------------------
# Ensure labels are factors with consistent levels
raw_data$true_labels <- factor(raw_data$true_labels, levels = c(0, 1))
raw_data$pred_labels <- factor(raw_data$pred_labels, levels = c(0, 1))

datasets <- unique(raw_data$Dataset)

# 4. Define Color Palettes
# ------------------------------------------------------------------------------
color_palettes <- list()
color_palettes[["Dataset2---Internal Testing"]] <- c("#FEFDFB", "#C1DEE6", "#79B5D0", "#1C93B5")
color_palettes[["Dataset3---Multicenter"]]       <- c("#FEFDFB", "#C5D5B8", "#97B786", "#3C7526")
color_palettes[["Dataset4---TCGA cohort"]]       <- c("#FEFDFB", "#FFE2BC", "#FFBE67", "#F2A312")

# 5. Plotting Function
# ------------------------------------------------------------------------------
plot_cm_square <- function(cm_table, dataset_name, palette_colors) {
  
  # Calculate percentages for heatmap intensity
  total_count <- sum(cm_table$Freq)
  cm_table$Freq_percentage <- cm_table$Freq / total_count
  
  # Generate Heatmap
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    # Tile layer with standardized linewidth
    geom_tile(color = "white", linewidth = 1.2) + 
    
    # Add count and percentage labels
    geom_text(aes(label = paste0(Freq, "\n", "(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 5) +
    
    scale_fill_gradientn(colors = palette_colors,
                         limits = c(0, 1),
                         labels = scales::percent,
                         name = NULL) +
    
    scale_x_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    scale_y_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    
    labs(title = dataset_name, x = "Predicted", y = "Actual") +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11, color = "black"),
      panel.grid = element_blank(), 
      legend.position = "right"
    ) +
    
    coord_fixed()
  
  return(p)
}

# 6. Generate Plots
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # Calculate Confusion Matrix using caret
  cm_obj <- confusionMatrix(sub_data$pred_labels, sub_data$true_labels)
  cm_df <- as.data.frame(cm_obj$table)
  
  # Assign palette (Fallback to grayscale if missing)
  current_palette <- color_palettes[[ds]]
  if (is.null(current_palette)) current_palette <- c("#FFFFFF", "#CCCCCC", "#666666", "#000000")
  
  plot_list[[ds]] <- plot_cm_square(cm_df, ds, current_palette)
}

# 7. Export Results
# ------------------------------------------------------------------------------
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(output_dir, "Figure5_Confusion_Matrices.pdf"), 
       plot = combined_plot, width = 15, height = 5, device = "pdf")

message("Analysis complete. Output saved to /results/ directory.")


# ==============================================================================
