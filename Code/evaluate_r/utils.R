# ==============================================================================
# Script Name: utils.R
# Description: Helper functions and environment setup for evaluation scripts.
# ==============================================================================

setup_environment <- function() {
  # Automatically install and load required packages
  required_packages <- c("ggplot2", "pROC", "dplyr", "gridExtra", "here")
  new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  
  library(ggplot2)
  library(pROC)
  library(dplyr)
  library(gridExtra)
  library(here)
  
  message("Environment setup complete. Packages loaded.")
}
