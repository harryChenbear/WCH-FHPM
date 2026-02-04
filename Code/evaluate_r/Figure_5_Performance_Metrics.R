# Script Name: Figure_5_Performance_Metrics.R
# Description: Calculate diagnostic metrics (Sensitivity, Specificity, PPV, NPV)
#              with 95% Confidence Intervals and generate a forest plot.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(scales)
library(here)

# Load data using relative path
data_path <- here("data", "Figure 5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

df <- read.csv(data_path)

# 2. Metric Calculation
# ------------------------------------------------------------------------------
# Confusion matrix components
TP <- sum(df$true_labels == 1 & df$pred_labels == 1)
TN <- sum(df$true_labels == 0 & df$pred_labels == 0)
FP <- sum(df$true_labels == 0 & df$pred_labels == 1)
FN <- sum(df$true_labels == 1 & df$pred_labels == 0)

# Function for Exact Binomial Confidence Intervals
get_ci <- function(success, total, metric_name) {
  test_result <- binom.test(success, total)
  data.frame(
    Metric = metric_name,
    Value = test_result$estimate[[1]],
    Lower = test_result$conf.int[1],
    Upper = test_result$conf.int[2]
  )
}

# Compile data frame
plot_data <- rbind(
  get_ci(TP, TP + FN, "Sensitivity"),
  get_ci(TN, TN + FP, "Specificity"),
  get_ci(TP, TP + FP, "PPV"),
  get_ci(TN, TN + FN, "NPV")
)

# Set factor levels for plotting order
plot_data$Metric <- factor(plot_data$Metric, 
                           levels = c("NPV", "PPV", "Specificity", "Sensitivity"))

# 3. Visualization
# ------------------------------------------------------------------------------
custom_color <- "#900C3F"

p <- ggplot(plot_data, aes(x = Value, y = Metric)) +
  # Confidence interval bars
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, 
                 color = custom_color, linewidth = 1) +
  
  # Point estimates
  geom_point(size = 5, color = custom_color) +
  
  # Percentage labels
  geom_text(aes(label = sprintf("%.1f%%", Value * 100)), 
            vjust = -1.5, color = custom_color, size = 4.5) +
  
  # Axis formatting
  scale_x_continuous(labels = percent_format(accuracy = 1), 
                     limits = c(0, 1.1), 
                     breaks = seq(0, 1, 0.2)) +
  
  labs(x = "Value (95% CI)", y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

# 4. Final Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(file.path(results_dir, "Figure_5_ForestPlot.pdf"), 
       plot = p, width = 8, height = 6)

message("Analysis complete. Forest plot saved to /results/ directory.")


# ==============================================================================
