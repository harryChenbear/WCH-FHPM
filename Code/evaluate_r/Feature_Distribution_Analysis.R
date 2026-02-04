# Script Name: Feature_Distribution_Analysis.R
# Description: Generates Violin and Boxplots with Wilcoxon statistical tests.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(here)

# Load data using relative path
data_path <- here("data", "Figure 4C.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 4C.csv' is in the /data/ directory.")
}

demo_data <- read.table(data_path, sep = ",", header = TRUE)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
# Reshape data for faceted plotting
plot_data <- demo_data %>%
  pivot_longer(cols = -1, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = names(demo_data)[-1]))

group_col <- names(demo_data)[1]

# 3. Statistical Analysis
# ------------------------------------------------------------------------------
stat_test <- plot_data %>%
  group_by(variable) %>%
  summarise(p_value = wilcox.test(value ~ .data[[group_col]])$p.value) %>%
  mutate(
    stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value >= 0.05 ~ "NS",
      TRUE ~ ""
    ),
    # Set significance label position 10% above max value
    y_pos = sapply(variable, function(var) {
      max_val <- max(plot_data$value[plot_data$variable == var], na.rm = TRUE)
      max_val * 1.1 
    })
  )

# 4. Visualization
# ------------------------------------------------------------------------------
p <- ggplot(plot_data, aes(x = .data[[group_col]], y = value, fill = .data[[group_col]])) +
  # Distribution layers
  geom_violin(alpha = 0.6, width = 0.8, trim = TRUE) + 
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, color = "black") + 
  
  # Significance annotation
  geom_text(data = stat_test, 
            aes(x = 1.5, y = y_pos, label = stars), 
            inherit.aes = FALSE, size = 5) +
  
  # Faceting: 5 columns for organized layout
  facet_wrap(~variable, scales = "free_y", ncol = 5, 
             labeller = labeller(variable = function(x) str_replace_all(x, "\\.", " "))) + 
  
  scale_fill_brewer(palette = "Set1") + 
  
  # Professional theme customization
  theme_minimal(base_size = 14) + 
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    axis.ticks = element_line(color = "black")
  ) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + 
  labs(x = "Group", y = "Feature Value", title = "Feature Distributions by Group") +
  coord_cartesian(ylim = c(0, NA))

# 5. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(filename = file.path(results_dir, "Figure_4C_Distributions.pdf"), 
       plot = p, width = 14, height = 8)

message("Analysis complete. Plot saved to /results/ directory.")


# ==============================================================================
