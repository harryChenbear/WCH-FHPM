# Script Name: Figure_3_Combined_Analysis.R
# Description: ROC curve analysis with overlay of pathologist performance metrics.
# Author: Jinge Zhao
# ==============================================================================

# 1. Setup Environment & Data
# ------------------------------------------------------------------------------
set.seed(123) 

# Load necessary libraries
required_packages <- c("pROC", "ggplot2", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(pROC)
library(here)

# Load dataset using relative path
data_path <- here("data", "Figure 3.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 3.csv' is in the /data/ directory.")
}

demo_data <- read.table(data_path, sep=",", header = TRUE)

# Identify the correct column for Fusion Strategy
strategy_col <- grep("Fusion[._ ]Strategy", colnames(demo_data), value = TRUE)
if (length(strategy_col) == 0) {
  stop("Error: Fusion Strategy column not found.")
}

unique_strategies <- unique(demo_data[[strategy_col]])
num_strategies <- length(unique_strategies)

# 2. Initialize Output Device
# ------------------------------------------------------------------------------
output_path <- here("results", "Figure3_Combined.pdf")
if (!dir.exists(here("results"))) dir.create(here("results"))

# Initialize PDF with dynamic height based on the number of strategies
pdf(output_path, width = 8, height = 6 * num_strategies)

# Set global plotting parameters: Square aspect ratio and margins
par(mfrow = c(num_strategies, 1), pty = "s", mar = c(5, 5, 4, 2))

# 3. Plotting Loop
# ------------------------------------------------------------------------------
for (strat in unique_strategies) {
  
  # Subset data for the current strategy
  current_data <- demo_data[demo_data[[strategy_col]] == strat, ]
  
  if(nrow(current_data) < 2) next
  
  # Calculate ROC model
  roc_model <- roc(current_data$true_labels, current_data$pred, quiet = TRUE)
  
  # --- Plot ROC Base ---
  plot(roc_model, 
       col = "#98d4da", 
       lwd = 4, 
       xlab = "Specificity", 
       ylab = "Sensitivity", 
       main = paste("Strategy:", strat), 
       cex.main = 1.8, 
       cex.lab = 1.5, 
       cex.axis = 1.2, 
       font.main = 2, 
       font.lab = 2,
       legacy.axes = TRUE,
       xlim = c(1, 0), 
       ylim = c(0, 1))
  
  grid(col = "gray85", lty = "dotted", lwd = 1)
  
  # --- Overlay Pathologist Performance ---
  for (i in 1:7) {
    # Extract binary predictions
    p_first <- current_data[[paste0("pathologist", i, "_first")]]
    p_second <- current_data[[paste0("pathologist", i, "_second")]]
    
    # Internal function to calculate Specificity and Sensitivity
    calc_perf <- function(pred, truth) {
      tp <- sum(pred == 1 & truth == 1); fp <- sum(pred == 1 & truth == 0)
      fn <- sum(pred == 0 & truth == 1); tn <- sum(pred == 0 & truth == 0)
      sens <- if((tp+fn)>0) tp/(tp+fn) else 0
      spec <- if((tn+fp)>0) tn/(tn+fp) else 0
      return(c(spec, sens))
    }
    
    res_first <- calc_perf(p_first, current_data$true_labels)
    res_second <- calc_perf(p_second, current_data$true_labels)
    
    # Color mapping for Junior (1-4) and Senior (5-7) groups
    if (i <= 4) {
      col_1 <- "#377eb8"; col_2 <- "#7a5195"; arr_col <- "#377eb8"
    } else {
      col_1 <- "#f5c656"; col_2 <- "#ed6a6a"; arr_col <- "#f5c656"
    }
    
    # Plot points and directional arrows
    points(res_first[1], res_first[2], col = col_1, pch = 19, cex = 2)
    points(res_second[1], res_second[2], col = col_2, pch = 17, cex = 2)
    arrows(res_first[1], res_first[2], res_second[1], res_second[2], 
           col = arr_col, length = 0.12, angle = 25, lwd = 2)
  }
  
  # --- Add Legend ---
  legend("bottomright", 
         legend = c("WCH-FHPM", 
                    "Junior pathologists (before)", 
                    "Junior pathologists (after)", 
                    "Senior pathologists (before)", 
                    "Senior pathologists (after)"), 
         col = c("#98d4da", "#377eb8", "#7a5195", "#f5c656", "#ed6a6a"), 
         lwd = 3, 
         pch = c(NA, 19, 17, 19, 17), 
         pt.cex = 1.5, 
         bty = "n", 
         text.font = 2, 
         cex = 1.0, 
         xpd = TRUE, 
         inset = c(0.02, 0.02))
}

# 4. Finalize
# ------------------------------------------------------------------------------
dev.off() 
message("Analysis complete. Combined PDF saved to 'results' directory.")


# ==============================================================================
