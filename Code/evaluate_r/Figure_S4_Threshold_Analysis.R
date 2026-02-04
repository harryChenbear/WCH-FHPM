# Script Name: Figure_S4_Threshold_Analysis.R
# Description: Multi-dataset confusion matrix analysis across varying cutoffs.
#              - Rows: Datasets; Columns: Probability Cutoffs (0.3 - 0.8).
#              - Independent color scales and legends per dataset.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(caret)
library(gridExtra)
library(grid)
library(scales) 
library(cowplot)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S4.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

data <- read.csv(data_path)
data$true_labels <- as.numeric(as.character(data$true_labels))
data$pred_probs  <- as.numeric(as.character(data$pred_probs))

datasets <- unique(data$Dataset)
cutoffs <- seq(0.3, 0.8, by = 0.1)

# 2. Color Palette Configuration
# ------------------------------------------------------------------------------
color_palettes <- list(
  "Dataset1---Training"         = c("#FEFDFB","#D4D1E6","#A29CCC","#7673A6"),
  "Dataset2---Internal Testing" = c("#FEFDFB","#C1DEE6","#79B5D0","#1C93B5"),
  "Dataset3---Multicenter"      = c("#FEFDFB","#C5D5B8","#97B786","#3C7526"),
  "Dataset4---TCGA cohort"      = c("#FEFDFB","#FFE2BC","#FFBE67","#F2A312")
)

# 3. Plotting Functions
# ------------------------------------------------------------------------------
# Function to generate individual confusion matrix tile
plot_confusion_matrix <- function(cm, cutoff, palette_colors, is_first_col, is_first_row) {
  cm_table <- as.data.frame(cm$table)
  cm_table$Freq_percentage <- cm_table$Freq / sum(cm_table$Freq)
  
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    geom_tile(color = "white", linewidth = 0.5) + 
    geom_text(aes(label = paste0(Freq, "\n(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 3) + 
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), name = NULL) +
    scale_x_discrete(labels = c("0" = "non-FH-dRCC", "1" = "FH-dRCC")) + 
    scale_y_discrete(labels = c("0" = "non-FH-dRCC", "1" = "FH-dRCC")) + 
    labs(
      title = if(is_first_row) paste0("Cutoff: ", cutoff) else NULL, 
      x = "Predicted", 
      y = if(is_first_col) "Actual" else NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = if(is_first_col) element_text(size = 9, face = "bold") else element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ) +
    coord_fixed()
  
  return(p)
}

# Function to extract colorbar legend for each row
get_row_legend <- function(palette_colors) {
  p <- ggplot(data.frame(x=1, y=1, val=0.5), aes(x, y, fill=val)) +
    geom_tile() +
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), 
                         labels = scales::percent, name = "Frequency") +
    guides(fill = guide_colorbar(barheight = unit(3.5, "cm"), barwidth = unit(0.4, "cm"),
                                 title.position = "top", title.hjust = 0.5)) +
    theme_minimal() +
    theme(legend.position = "right", 
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))
  return(get_legend(p))
}

# 4. Processing Loop
# ------------------------------------------------------------------------------
plot_list <- list()
legend_list <- list() 
row_idx <- 0

for (ds in datasets) {
  row_idx <- row_idx + 1
  sub_data <- subset(data, Dataset == ds)
  current_palette <- color_palettes[[ds]] %||% c("#FFFFFF", "#CCCCCC", "#000000")
  
  legend_list[[ds]] <- get_row_legend(current_palette)
  
  col_idx <- 0
  for (cutoff in cutoffs) {
    col_idx <- col_idx + 1
    pred_factor <- factor(ifelse(sub_data$pred_probs >= cutoff, 1, 0), levels = c(0, 1))
    true_factor <- factor(sub_data$true_labels, levels = c(0, 1))
    
    cm <- confusionMatrix(pred_factor, true_factor)
    p <- plot_confusion_matrix(cm, cutoff, current_palette, (col_idx == 1), (row_idx == 1))
    plot_list[[paste0(ds, "_", cutoff)]] <- p
  }
}

# 5. Assembly and Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Construct sub-grids
main_grid <- arrangeGrob(grobs = plot_list, nrow = length(datasets), ncol = length(cutoffs))
legend_col <- arrangeGrob(grobs = legend_list, ncol = 1)
row_labels <- arrangeGrob(grobs = lapply(datasets, function(name) {
  textGrob(strsplit(as.character(name), "---")[[1]][1], rot = -90, gp = gpar(fontsize = 12, fontface = "bold"))
}), ncol = 1)

# Final Arrangement
final_plot <- grid.arrange(main_grid, legend_col, row_labels, ncol = 3, widths = c(10, 2, 0.5))

ggsave(filename = file.path(results_dir, "Figure_S4_Threshold_Analysis.pdf"), 
       plot = final_plot, width = length(cutoffs) * 2.2 + 2.5, height = length(datasets) * 2.4, device = "pdf")

message("Analysis complete. Results saved to /results/.")


# ==============================================================================
