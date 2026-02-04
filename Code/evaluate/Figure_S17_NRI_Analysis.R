# Script Name: Figure_S17_NRI_Analysis.R
# Description: Analysis of Net Reclassification Index (NRI) and Confusion Matrices
#              to quantify diagnostic improvement with WCHFHPM model support.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(here)

# Use relative path for data import
data_path <- here("data", "Figure S17.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

demo_data <- read.table(data_path, sep = ",", header = TRUE)
true_labels <- demo_data$true_labels

# Initialize list to store plot objects for all readers
all_plots <- list()

# 2. NRI Calculation and Plotting Loop (Readers 1 to 7)
# ------------------------------------------------------------------------------
for (i in 1:7) {
  
  message(paste("Analyzing performance for Reader", i, "..."))
  
  # Extract prediction vectors
  pred_without <- demo_data[[paste0("Human", i, "_Without_WCHFHPM")]]
  pred_with    <- demo_data[[paste0("Human", i, "_With_WCHFHPM")]]
  
  # --- Part A: Confusion Matrices (Comparison within Ground Truth subsets) ---
  
  # Subset 1: True FHdRCC cases
  idx_1 <- which(true_labels == 1)
  cm_1  <- table(First = factor(pred_without[idx_1], levels=0:1), 
                 Second = factor(pred_with[idx_1], levels=0:1))
  
  # Subset 2: True non-FHdRCC cases
  idx_0 <- which(true_labels == 0)
  cm_0  <- table(First = factor(pred_without[idx_0], levels=0:1), 
                 Second = factor(pred_with[idx_0], levels=0:1))
  
  # --- Part B: NRI Calculation ---
  
  # Reclassification logic
  improved_tp <- sum(pred_without == 0 & pred_with == 1 & true_labels == 1) # FN -> TP
  worsened_fn <- sum(pred_without == 1 & pred_with == 0 & true_labels == 1) # TP -> FN
  
  improved_tn <- sum(pred_without == 1 & pred_with == 0 & true_labels == 0) # FP -> TN
  worsened_fp <- sum(pred_without == 0 & pred_with == 1 & true_labels == 0) # TN -> FP
  
  # Probabilities of reclassification
  p_event_up    <- (improved_tp - worsened_fn) / sum(true_labels == 1)
  p_nonevent_up <- (improved_tn - worsened_fp) / sum(true_labels == 0)
  
  total_nri <- p_event_up + p_nonevent_up
  
  # --- Part C: Visualization ---
  
  # Prepare data for Heatmaps
  cm1_df <- as.data.frame(as.table(cm_1))
  levels(cm1_df$First) <- levels(cm1_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  cm0_df <- as.data.frame(as.table(cm_0))
  levels(cm0_df$First) <- levels(cm0_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  # Plot 1: FHdRCC Confusion Matrix
  p1 <- ggplot(cm1_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    labs(title = "GT: FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 2: non-FHdRCC Confusion Matrix
  p2 <- ggplot(cm0_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15") +
    labs(title = "GT: non-FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 3: NRI Component Analysis
  nri_df <- data.frame(
    Category = c("FHdRCC Imp.", "non-FHdRCC Imp.", "Total NRI"),
    Value = c(p_event_up, p_nonevent_up, total_nri)
  )
  
  p3 <- ggplot(nri_df, aes(x = Category, y = Value, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    labs(title = "NRI Results", x = NULL, y = "Value") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine reader-specific charts
  all_plots[[i]] <- p1 + p2 + p3 + plot_layout(ncol = 3)
}

# 3. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Vertically stack all reader panels
final_fig <- wrap_plots(all_plots, ncol = 1)

ggsave(filename = file.path(results_dir, "Figure_S17_NRI_Analysis.pdf"), 
       plot = final_fig, width = 12, height = 4 * 7, limitsize = FALSE)

message("Analysis complete. Combined PDF exported to /results/.")


# ==============================================================================
