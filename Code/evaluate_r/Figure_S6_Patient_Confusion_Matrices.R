# Script Name: Figure_S6_Patient_Confusion_Matrices.R
# Description: Generate Patient-Level Square Confusion Matrices for Figure S6.
#              Aggregation Rule: Max pred_labels per Patient ID.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(caret)
library(gridExtra)
library(scales)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing (Aggregation)
# ------------------------------------------------------------------------------
# Identify Patient ID column (matching "Pts")
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (!is.na(pts_col)) colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Ensure numeric types
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
if(!"pred_labels" %in% colnames(raw_data)) stop("Column 'pred_labels' not found.")
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# Aggregate: Take MAX pred_label for each patient (Patient-level diagnosis)
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>%
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to factors for Confusion Matrix
patient_data$true_labels <- factor(patient_data$true_labels, levels = c(0, 1))
patient_data$pred_labels <- factor(patient_data$pred_labels, levels = c(0, 1))

datasets <- unique(patient_data$Dataset)

# 3. Define Color Palettes
# ------------------------------------------------------------------------------
color_palettes <- list(
  "Dataset2---Internal Testing" = c("#FEFDFB", "#C1DEE6", "#79B5D0", "#1C93B5"),
  "Dataset3---Multicenter"      = c("#FEFDFB", "#C5D5B8", "#97B786", "#3C7526"),
  "Dataset4---TCGA cohort"      = c("#FEFDFB", "#FFE2BC", "#FFBE67", "#F2A312")
)

# 4. Plotting Function
# ------------------------------------------------------------------------------
plot_cm_square <- function(cm_table, dataset_name, palette_colors) {
  
  # Calculate percentages
  total_count <- sum(cm_table$Freq)
  cm_table$Freq_percentage <- cm_table$Freq / total_count
  
  # Generate Heatmap
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    geom_tile(color = "white", linewidth = 1.2) + 
    geom_text(aes(label = paste0(Freq, "\n(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 5) +
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), labels = scales::percent, name = NULL) +
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

# 5. Generation and Export
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(patient_data, Dataset == ds)
  cm_obj <- confusionMatrix(sub_data$pred_labels, sub_data$true_labels)
  cm_df <- as.data.frame(cm_obj$table)
  
  current_palette <- color_palettes[[ds]] %||% c("#FFFFFF", "#CCCCCC", "#000000")
  plot_list[[ds]] <- plot_cm_square(cm_df, ds, current_palette)
}

# Arrange and save
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure_S6_Patient_Confusion_Matrices.pdf"), 
       plot = combined_plot, width = 15, height = 5, device = "pdf")

message("Patient-level matrices saved to /results/.")


# ==============================================================================
