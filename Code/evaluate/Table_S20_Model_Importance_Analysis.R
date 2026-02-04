# Script Name: Table_S20_Model_Importance_Analysis.R
# Description: Train logistic regression models, extract statistical equations, 
#              and visualize standardized feature importance.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(caret)      # For data partitioning and pre-processing
library(ggplot2)    # For visualizations
library(gridExtra)  # For plot arrangement
library(speedglm)   # For efficient model fitting
library(here)       # For relative path management

# Load dataset using relative path
data_path <- here("data", "Table S20.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

data <- read.table(data_path, sep = ",", header = TRUE)

# Standardize column names: outcome (label), AI prediction (pred), clinical features (f1-f5)
colnames(data) <- c("label", "pred", "f1", "f2", "f3", "f4", "f5")

# 2. Helper Functions
# ------------------------------------------------------------------------------

# Z-score standardization based on training set statistics
standardize <- function(train_data, test_data, features) {
  scaler <- preProcess(train_data[, features], method = c("center", "scale"))
  
  train_scaled <- predict(scaler, train_data[, features])
  test_scaled <- predict(scaler, test_data[, features])
  
  list(
    train = cbind(train_scaled, label = train_data$label),
    test = cbind(test_scaled, label = test_data$label)
  )
}

# Extract model equation and p-values
get_logistic_info <- function(model) {
  coef_values <- coef(model)
  vcov_matrix <- vcov(model)
  std_err     <- sqrt(diag(vcov_matrix))
  z_values    <- coef_values / std_err
  p_values    <- 2 * pnorm(-abs(z_values))
  
  formula_terms <- names(coef_values)[-1]
  
  list(
    formula = paste("logit(p) =", round(coef_values[1], 3), 
                    "+", paste(round(coef_values[-1], 3), formula_terms, collapse = " + ")),
    p_values = paste(sapply(seq_along(formula_terms), function(i) {
      paste(formula_terms[i], "p-value =", format.pval(p_values[i+1], digits = 3))
    }), collapse = "\n")
  )
}

# Generate importance plot using absolute coefficient values
plot_importance <- function(model, title) {
  coef_df <- data.frame(
    feature = names(coef(model))[-1],
    importance = abs(coef(model)[-1])
  )
  
  ggplot(coef_df, aes(x = reorder(feature, importance), y = importance)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    coord_flip() +
    labs(title = title, x = "Features", y = "Importance (|Coefficient|)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# 3. Model Training Pipeline
# ------------------------------------------------------------------------------
set.seed(123) 
trainIndex <- createDataPartition(data$label, p = 0.7, list = FALSE)
train_raw  <- data[trainIndex, ]
test_raw   <- data[-trainIndex, ]

# Apply Z-score standardization (essential for ranking feature importance)
feature_list <- c("f1", "f2", "f3", "f4", "f5")
scaled_data <- standardize(train_raw, test_raw, feature_list)
train_full  <- scaled_data$train
train_reduced <- train_full[, c("f1", "f2", "label")]

# Fit models
# Model 1: Reduced set
model_reduced <- speedglm(label ~ ., data = train_reduced, family = binomial(), method = "Cholesky")

# Model 2: Full set
model_full <- speedglm(label ~ ., data = train_full, family = binomial(), method = "Cholesky")

# 4. Results Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Print statistics to console
mod1_info <- get_logistic_info(model_reduced)
mod2_info <- get_logistic_info(model_full)

cat("\n--- Reduced Model ---\n", mod1_info$formula, "\n", mod1_info$p_values, "\n")
cat("\n--- Full Model ---\n", mod2_info$formula, "\n", mod2_info$p_values, "\n")

# Save visualizations to PDF
output_path <- file.path(results_dir, "Table_S20_Feature_Importance.pdf")
pdf(file = output_path, width = 10, height = 5)

grid.arrange(
  plot_importance(model_reduced, "Importance: Reduced Model"),
  plot_importance(model_full, "Importance: Full Model"),
  ncol = 2
)

dev.off()
message("Feature importance analysis complete. Output saved to: ", results_dir)