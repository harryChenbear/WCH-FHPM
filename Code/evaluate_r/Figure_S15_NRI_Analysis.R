# Script Name: Figure_S15_NRI_Analysis.R
# Description: Analysis of Net Reclassification Index (NRI) and Confusion Matrices
#              comparing reader performance "Without" vs "With" WCHFHPM model.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(here)

# Load Data using relative path
data_path <- here("data", "Figure S15.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

demo_data <- read.table(data_path, sep = ",", header = TRUE)
true_labels <- demo_data$true_labels

# Initialize storage for diagnostic plots
all_plots <- list()

# 2. NRI Calculation and Plotting Loop (Readers 1-7)
# ------------------------------------------------------------------------------
for (i in 1:7) {
  
  message(paste("Calculating NRI for Reader", i, "..."))
  
  # Extract specific reader predictions
  without_ai <- demo_data[[paste0("Human", i, "_Without_WCHFHPM")]]
  with_ai    <- demo_data[[paste0("Human", i, "_With_WCHFHPM")]]
  
  # --- Part A: Confusion Matrices for Ground Truth Subsets ---
  
  # 1. Confusion matrix for FHdRCC cases (True Positive subset)
  idx_1 <- which(true_labels == 1)
  cm_1  <- table(First = factor(without_ai[idx_1], levels=0:1), 
                 Second = factor(with_ai[idx_1], levels=0:1))
  
  # 2. Confusion matrix for non-FHdRCC cases (True Negative subset)
  idx_0 <- which(true_labels == 0)
  cm_0  <- table(First = factor(without_ai[idx_0], levels=0:1), 
                 Second = factor(with_ai[idx_0], levels=0:1))
  
  # --- Part B: Net Reclassification Index (NRI) Logic ---
  
  # Reclassification counts
  improved_tp <- sum(without_ai == 0 & with_ai == 1 & true_labels == 1) # FN -> TP
  worsened_fn <- sum(without_ai == 1 & with_ai == 0 & true_labels == 1) # TP -> FN
  
  improved_tn <- sum(without_ai == 1 & with_ai == 0 & true_labels == 0) # FP -> TN
  worsened_fp <- sum(without_ai == 0 & with_ai == 1 & true_labels == 0) # TN -> FP
  
  # Probability of improvement
  p_event_up    <- (improved_tp - worsened_fn) / sum(true_labels == 1)
  p_nonevent_up <- (improved_tn - worsened_fp) / sum(true_labels == 0)
  
  total_nri <- p_event_up + p_nonevent_up
  
  # --- Part C: Visualization ---
  
  # Data preparation for heatmaps
  cm1_df <- as.data.frame(as.table(cm_1))
  levels(cm1_df$First) <- levels(cm1_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  cm0_df <- as.data.frame(as.table(cm_0))
  levels(cm0_df$First) <- levels(cm0_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  # Heatmap 1: FHdRCC Cases
  p1 <- ggplot(cm1_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    labs(title = "Ground truth: FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Heatmap 2: non-FHdRCC Cases
  p2 <- ggplot(cm0_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15") +
    labs(title = "Ground truth: non-FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Bar Chart: NRI Components
  nri_df <- data.frame(
    Category = c("FHdRCC Improvement", "non-FHdRCC Improvement", "Total NRI"),
    Value = c(p_event_up, p_nonevent_up, total_nri)
  )
  
  p3 <- ggplot(nri_df, aes(x = Category, y = Value, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    labs(title = "NRI Components", x = NULL, y = "Index Value") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine reader-specific plots
  all_plots[[i]] <- p1 + p2 + p3 + plot_layout(ncol = 3)
}

# 3. Final Assembly and Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

final_plot <- wrap_plots(all_plots, ncol = 1) # Vertical stack for multiple readers

# Save as high-resolution PDF
ggsave(filename = file.path(results_dir, "Figure_S15_NRI_Analysis.pdf"), 
       plot = final_plot, width = 12, height = 4 * 7, limitsize = FALSE)

message("NRI Analysis complete. Results saved to /results/.")


# ==============================================================================
