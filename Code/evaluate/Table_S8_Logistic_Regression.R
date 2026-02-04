# Script Name: Table_S8_Logistic_Regression.R
# Description: Performs univariate logistic regression for multiple predictors.
#              Calculates Odds Ratios (OR), 95% CI, and P-values.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(here)

# Load data using relative path
data_path <- here("data", "Table S8.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

# 'check.names = FALSE' preserves spaces in clinical feature names
data <- read.csv(data_path, check.names = FALSE)

# 2. Variable Definition
# ------------------------------------------------------------------------------
# Predictors: First 11 columns; Target: 12th column
predictor_cols <- colnames(data)[1:11]
target_col <- colnames(data)[12]

# Initialize storage for results
results <- data.frame(
  Feature = character(),
  OR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# 3. Univariate Logistic Regression Loop
# ------------------------------------------------------------------------------
message("Starting univariate logistic regression...")

for (var in predictor_cols) {
  
  # Construct formula using backticks to handle special characters/spaces
  formula_str <- paste0("`", target_col, "` ~ `", var, "`")
  
  # Fit Generalized Linear Model (Binomial family)
  model <- glm(as.formula(formula_str), data = data, family = binomial)
  model_sum <- summary(model)
  
  # Check for coefficient validity (avoid constant variable errors)
  if (nrow(model_sum$coefficients) > 1) {
    beta <- coef(model)[2]
    p_val <- model_sum$coefficients[2, 4]
    
    # Calculate Odds Ratio
    OR_val <- exp(beta)
    
    # Calculate 95% Confidence Intervals
    # Defaults to Wald CI if Profile Likelihood fails due to convergence issues
    ci <- tryCatch({
      exp(confint(model)[2, ]) 
    }, error = function(e) {
      exp(confint.default(model)[2, ]) 
    })
    
    ci_low <- ci[1]
    ci_high <- ci[2]
  } else {
    OR_val <- NA
    ci_low <- NA
    ci_high <- NA
    p_val  <- NA
  }
  
  # Append to results
  results <- rbind(results, data.frame(
    Feature = var,
    OR = OR_val,
    CI_Lower = ci_low,
    CI_Upper = ci_high,
    P_Value = p_val
  ))
}

# 4. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(results, file.path(results_dir, "Table_S8_Logistic_Results.csv"), row.names = FALSE)

print(results)
message("Analysis complete. Results saved to /results/.")


# ==============================================================================
