# Script Name: Figure_4D_Model_Comparison.R
# Description: Comparative ROC analysis between AI (WCH-FHPM) and 
#              Morphological Models.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(ggplot2)
library(here)

# Load data using relative path
data_path <- here("data", "Figure 4D.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

df <- read.csv(data_path)

# 2. ROC Curve and AUC Calculation
# ------------------------------------------------------------------------------
# Calculate ROC and AUC for each model
# Note: Ensure the column name 'true_lables' matches your CSV header
roc_ai <- roc(df$true_lables, df$pred_probs, quiet = TRUE)
auc_ai <- auc(roc_ai)

roc_m1 <- roc(df$true_lables, df$Model1_probs, quiet = TRUE)
auc_m1 <- auc(roc_m1)

roc_m2 <- roc(df$true_lables, df$Model2_probs, quiet = TRUE)
auc_m2 <- auc(roc_m2)

# Log AUC results
cat(sprintf("AUC values:\nAI: %.3f\nModel2: %.3f\nModel1: %.3f\n", auc_ai, auc_m2, auc_m1))

# 3. Data Preparation for Plotting
# ------------------------------------------------------------------------------
# Consolidate ROC data into a single data frame for ggplot
plot_data <- rbind(
  data.frame(spec = roc_ai$specificities, sens = roc_ai$sensitivities, 
             Model = sprintf("WCH-FHPM (AUC=%.3f)", auc_ai)),
  data.frame(spec = roc_m1$specificities, sens = roc_m1$sensitivities, 
             Model = sprintf("Morphological model-1 (AUC=%.3f)", auc_m1)),
  data.frame(spec = roc_m2$specificities, sens = roc_m2$sensitivities, 
             Model = sprintf("Morphological model-2 (AUC=%.3f)", auc_m2))
)

# Enforce factor levels for legend ordering
plot_data$Model <- factor(plot_data$Model, levels = c(
  sprintf("WCH-FHPM (AUC=%.3f)", auc_ai),
  sprintf("Morphological model-1 (AUC=%.3f)", auc_m1),
  sprintf("Morphological model-2 (AUC=%.3f)", auc_m2)
))

# 4. Visualization
# ------------------------------------------------------------------------------
custom_colors <- c("#993344", "#CCBBCC", "#4477AA") 

p <- ggplot(plot_data, aes(x = spec, y = sens, color = Model)) +
  # Diagonal reference line
  geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
               color = "lightgray", linetype = "solid", linewidth = 0.8) +
  
  # ROC Curves
  geom_path(linewidth = 1.2) +
  
  # Coordinate Configuration (Reversed Specificity)
  scale_x_reverse(name = "Specificity", limits = c(1.05, -0.05), expand = c(0,0)) +
  scale_y_continuous(name = "Sensitivity", limits = c(-0.05, 1.05), expand = c(0,0)) +
  
  scale_color_manual(values = custom_colors) +
  
  # Academic theme and layout customization
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.60, 0.25),
    legend.title = element_blank(),
    legend.background = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)
  )

# 5. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(file.path(results_dir, "Figure_4D_ROC_Comparison.pdf"), 
       plot = p, width = 6, height = 6)

message("Analysis complete. ROC plot saved to /results/ directory.")


# ==============================================================================
