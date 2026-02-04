# Script Name: Table_S19_Feature_Proportions.R
# Description: Generates symmetrical feature proportion plots comparing 
#              Stroma and Parenchyma characteristics across clinical groups.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 
library(here)

# Load data using relative path
data_path <- here("data", "Table S19.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

# check.names = FALSE ensures column names with spaces are preserved
demo_data <- read.csv(data_path, check.names = FALSE)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
# Define feature columns (assuming first 11 columns)
feature_cols <- names(demo_data)[1:11]
feature_order <- feature_cols

# Transform group labels for visualization
# Logic: true_lables == 2 (Case, shown upwards), true_lables == 1 (Control, shown downwards)
demo_data$group <- ifelse(demo_data$true_lables == 2, 1, 0)

# Retrieve unique comparison group names
comparison_groups <- unique(demo_data$`Comparison groups`)

# 3. Plotting Parameters
# ------------------------------------------------------------------------------
gap_height <- 0.15   # Vertical gap for central placeholders
box_width  <- 0.7    # Width of bars
stroma_end <- 4.5    # Position separator for tissue categories
num_features <- length(feature_order)
label_pos <- c(mean(1:4), mean(5:num_features))

# 4. Plotting Function
# ------------------------------------------------------------------------------
create_plot <- function(sub_data, group_title) {
  
  # Calculate Proportions
  feature_proportions <- sub_data %>%
    select(all_of(c("group", feature_cols))) %>%
    pivot_longer(cols = -group, names_to = "feature", values_to = "value") %>%
    group_by(group, feature) %>%
    summarise(proportion = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    mutate(
      feature = factor(feature, levels = feature_order),
      ymin = ifelse(group == 0, -(gap_height + proportion), gap_height),
      ymax = ifelse(group == 0, -gap_height, gap_height + proportion),
      label_y = ifelse(group == 0, -(gap_height + proportion) - 0.05, (gap_height + proportion) + 0.05)
    )
  
  # Statistical testing (t-test)
  p_values <- sub_data %>%
    select(all_of(c("group", feature_cols))) %>%
    pivot_longer(cols = -group, names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    summarise(p_value = t.test(value ~ group)$p.value, .groups = 'drop') %>%
    mutate(feature = factor(feature, levels = feature_order))
  
  # Merge data
  grouped_data <- feature_proportions %>% left_join(p_values, by = "feature")
  
  # Visualization
  p <- ggplot(grouped_data) +
    # Background shading for tissue categories
    annotate("rect", xmin = 0.5, xmax = stroma_end, ymin = -Inf, ymax = Inf, fill = "#F0F8FF", alpha = 0.3) +
    annotate("rect", xmin = stroma_end, xmax = num_features + 0.5, ymin = -Inf, ymax = Inf, fill = "#FFF0F5", alpha = 0.3) +
    
    # Central placeholders
    geom_rect(aes(xmin = as.numeric(feature) - box_width/2, 
                  xmax = as.numeric(feature) + box_width/2,
                  ymin = -gap_height, 
                  ymax = gap_height),
              fill = "white", color = "black", linetype = "dashed", linewidth = 0.5) +
    
    # Main Bars
    geom_rect(aes(xmin = as.numeric(feature) - box_width/2,
                  xmax = as.numeric(feature) + box_width/2,
                  ymin = ymin,
                  ymax = ymax,
                  fill = abs(ymax - ymin)), 
              color = NA) +
    
    # Category Headers
    annotate("segment", x = 0.5, xend = stroma_end, y = 1.15, yend = 1.15, color = "#1E90FF", linewidth = 1.2) +
    annotate("segment", x = stroma_end, xend = num_features + 0.5, y = 1.15, yend = 1.15, color = "#CD5C5C", linewidth = 1.2) +
    annotate("text", x = label_pos[1], y = 1.2, label = "Stroma characteristics", size = 5, fontface = "bold", color = "#1E90FF") +
    annotate("text", x = label_pos[2], y = 1.2, label = "Parenchyma characteristics", size = 5, fontface = "bold", color = "#CD5C5C") +
    
    # Scales and Themes
    scale_fill_gradient(low = "#B3E5FC", high = "#0288D1") +
    scale_x_discrete(limits = feature_order, labels = function(x) gsub("\\.", " ", x)) +
    scale_y_continuous(
      breaks = c(seq(-1, -gap_height, by = 0.2), seq(gap_height, 1, by = 0.2)),
      labels = function(y) {
        val <- abs(y) - gap_height
        val <- ifelse(val < 0, 0, val)
        paste0(round(val * 100), "%")
      },
      limits = c(-1.2, 1.3)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    ) +
    labs(title = group_title, x = NULL, y = "Proportion") +
    
    # Bar Labels
    geom_text(aes(x = feature, y = label_y, label = sprintf("%.1f%%", abs(proportion * 100))),
              size = 3.5, fontface = "bold", color = "black") +
    
    # Significance Stars
    geom_text(data = grouped_data %>% group_by(feature) %>% summarise(p_value = first(p_value)),
              aes(x = feature, y = 1.05,
                  label = case_when(p_value < 0.001 ~ '***',
                                    p_value < 0.01 ~ '**',
                                    p_value < 0.05 ~ '*',
                                    TRUE ~ 'ns')),
              inherit.aes = FALSE, size = 5, fontface = "bold", color = "black")
  
  return(p)
}

# 5. Execution
# ------------------------------------------------------------------------------
if(length(comparison_groups) == 0) stop("No comparison groups identified.")

plot_list <- lapply(comparison_groups, function(grp) {
  create_plot(demo_data %>% filter(`Comparison groups` == grp), grp)
})

final_plot <- wrap_plots(plot_list, ncol = 1)

# Export results
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(file.path(results_dir, "Table_S19_Feature_Proportions.png"), 
       final_plot, width = 12, height = 8 * length(comparison_groups), dpi = 300)

message("Analysis complete. Results exported to /results/.")


# ==============================================================================
