# Script Name: Figure_S9_Simulation_Analysis.R
# Description: Evaluates model stability via 1000 random sampling iterations.
#              Simulates a fixed prevalence (20 positive / 180 negative).
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)     
library(caret)    
library(knitr)    
library(ggplot2)  
library(reshape2) 
library(dplyr)    
library(here)     

# Load data using relative path
data_path <- here("data", "Figure S9.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

data <- read.table(data_path, sep = ",", header = TRUE) 

# 2. Simulation Configuration
# ------------------------------------------------------------------------------
set.seed(123) # For reproducibility
n_iterations <- 1000

# Pre-subset data by class to improve loop efficiency
data_1 <- subset(data, true_labels == 1)
data_0 <- subset(data, true_labels == 0)

# Initialize matrix to store results
results_mat <- matrix(NA, nrow = n_iterations, ncol = 6)
colnames(results_mat) <- c("Sensitivity", "Specificity", "PPV", "NPV", "AUROC", "Brier_Score")

# 3. Perform Simulation Loop
# ------------------------------------------------------------------------------
message("Starting 1000 iterations of random sampling simulation...")

for (i in 1:n_iterations) {
  # Randomly sample observations (20 positive, 180 negative)
  sample_1 <- data_1[sample(1:nrow(data_1), 20, replace = FALSE), ] 
  sample_0 <- data_0[sample(1:nrow(data_0), 180, replace = FALSE), ] 
  
  # Combine and shuffle
  sampled_data <- rbind(sample_1, sample_0)
  sampled_data <- sampled_data[sample(1:nrow(sampled_data)), ]
  
  y_true <- as.factor(sampled_data$true_labels)
  y_pred <- as.factor(sampled_data$pred_labels)
  
  # Confusion Matrix Metrics
  conf_matrix <- confusionMatrix(y_pred, y_true, positive = "1")
  
  sensitivity <- conf_matrix$byClass["Sensitivity"]
  specificity <- conf_matrix$byClass["Specificity"]
  ppv <- conf_matrix$byClass["Pos Pred Value"]
  npv <- conf_matrix$byClass["Neg Pred Value"]
  
  # AUROC calculation
  roc_obj <- roc(sampled_data$true_labels, sampled_data$pred_labels, quiet = TRUE)
  auroc <- as.numeric(auc(roc_obj))
  
  # Brier Score calculation
  brier_score <- mean((sampled_data$pred_labels - sampled_data$true_labels)^2)
  
  results_mat[i, ] <- c(sensitivity, specificity, ppv, npv, auroc, brier_score)
}

results <- as.data.frame(results_mat)

# 4. Summary and Export
# ------------------------------------------------------------------------------
results_summary <- results %>%
  summarise(across(everything(), list(Mean = ~mean(., na.rm = TRUE)))) %>%
  reshape2::melt() %>%
  rename(Metric = variable, Mean = value)

# Print Summary Table
print(kable(results_summary, format = "markdown"))

# 5. Visualization (Density Plot)
# ------------------------------------------------------------------------------
results_long <- reshape2::melt(results)
median_values <- results_long %>% 
  group_by(variable) %>% 
  summarise(value = median(value, na.rm = TRUE))

plot <- ggplot(results_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_vline(data = median_values, 
             aes(xintercept = value), color = "red", linetype = "dashed", linewidth = 1) +
  geom_label(data = median_values, aes(x = value, y = 0, label = sprintf("%.3f", value)), 
             color = "red", vjust = -0.5, size = 5, fontface = "bold", fill = "white", label.size = 0) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal(base_size = 15) +
  labs(title = "Distribution of Metrics over 1000 Iterations", x = "Value", y = "Density") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(linewidth = 0.5, linetype = "dashed", color = "grey80"),
    legend.position = "none"
  )

# Export Result
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(filename = file.path(results_dir, "Figure_S9_Simulation_Density.pdf"), 
       plot = plot, width = 12, height = 8, units = "in", dpi = 300)

message("Simulation complete. Density plot saved to /results/.")


# ==============================================================================
