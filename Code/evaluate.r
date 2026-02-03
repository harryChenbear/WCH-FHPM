# ==============================================================================
# Script Name: Figure_2_ROC_Analysis.R
# Description: Generate ROC curves for multiple datasets.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
# Automatically install and load required packages
required_packages <- c("ggplot2", "pROC", "dplyr", "gridExtra", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here) 

# 2. Data Import
# ------------------------------------------------------------------------------
# Load data from the 'data' directory
data_path <- here("data", "Figure 2.csv")

if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 2.csv' is in the /data/ directory.")
}

raw_data <- read.csv(data_path)

# 3. Data Preprocessing
# ------------------------------------------------------------------------------
datasets <- unique(raw_data$Dataset)
raw_data$Dataset <- factor(raw_data$Dataset, levels = datasets)

# Define color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 4. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list() 

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (AUC & CI) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  # Construct annotation string
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Coordinate Extraction ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  # --- Step C: Create Plot Object ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = custom_colors[ds], linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             size = 4, hjust = 0, lineheight = 1.2) +
    
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = ds) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[ds]] <- p
}

# 5. Export Results
# ------------------------------------------------------------------------------
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(output_dir, "Figure2_ROC_Analysis.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Output saved to: ", output_dir)




# ==============================================================================
# Script Name: Figure_2_Confusion_Matrices.R
# Description: Generate confusion matrices for datasets in Figure 2.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "caret", "gridExtra", "scales", "dplyr", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(caret)
library(gridExtra)
library(scales)
library(dplyr)
library(here)

# 2. Data Import
# ------------------------------------------------------------------------------
# Load data from the 'data' directory
data_path <- here("data", "Figure 2.csv")

if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 2.csv' is in the /data/ directory.")
}

raw_data <- read.csv(data_path)

# 3. Data Preprocessing
# ------------------------------------------------------------------------------
# Ensure labels are factors with consistent levels
raw_data$true_labels <- factor(raw_data$true_labels, levels = c(0, 1))
raw_data$pred_labels <- factor(raw_data$pred_labels, levels = c(0, 1))

datasets <- unique(raw_data$Dataset)

# 4. Define Color Palettes
# ------------------------------------------------------------------------------
color_palettes <- list()
color_palettes[["Dataset2---Internal Testing"]] <- c("#FEFDFB", "#C1DEE6", "#79B5D0", "#1C93B5")
color_palettes[["Dataset3---Multicenter"]]       <- c("#FEFDFB", "#C5D5B8", "#97B786", "#3C7526")
color_palettes[["Dataset4---TCGA cohort"]]       <- c("#FEFDFB", "#FFE2BC", "#FFBE67", "#F2A312")

# 5. Plotting Function
# ------------------------------------------------------------------------------
plot_cm_square <- function(cm_table, dataset_name, palette_colors) {
  
  # Calculate percentages for heatmap intensity
  total_count <- sum(cm_table$Freq)
  cm_table$Freq_percentage <- cm_table$Freq / total_count
  
  # Generate Heatmap
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    geom_tile(color = "white", linewidth = 1.2) + 
    
    # Add count and percentage labels
    geom_text(aes(label = paste0(Freq, "\n", "(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 5) +
    
    scale_fill_gradientn(colors = palette_colors,
                         limits = c(0, 1),
                         labels = scales::percent,
                         name = NULL) +
    
    scale_x_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    scale_y_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    
    labs(title = dataset_name, x = "Predicted", y = "Actual") +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11, color = "black"),
      panel.grid = element_blank(), 
      legend.position = "right"
    ) +
    
    coord_fixed()
  
  return(p)
}

# 6. Generate Plots
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # Calculate Confusion Matrix using caret
  cm_obj <- confusionMatrix(sub_data$pred_labels, sub_data$true_labels)
  cm_df <- as.data.frame(cm_obj$table)
  
  # Assign palette
  current_palette <- color_palettes[[ds]]
  if (is.null(current_palette)) current_palette <- c("#FFFFFF", "#CCCCCC", "#666666", "#000000")
  
  plot_list[[ds]] <- plot_cm_square(cm_df, ds, current_palette)
}

# 7. Export Results
# ------------------------------------------------------------------------------
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(output_dir, "Figure2_Confusion_Matrices.pdf"), 
       plot = combined_plot, width = 15, height = 5, device = "pdf")

message("Analysis complete. Output saved to: ", output_dir)


# ==============================================================================
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
# Script Name: Figure_4B_Feature_Analysis.R
# Description: Generates a symmetrical feature distribution plot with 
#              locked aspect ratios for center box visualization.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup & Data Loading
# ------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(grid)
library(here)

set.seed(123)

# Load data using relative path
data_path <- here("data", "Figure 4B.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

data <- read.csv(data_path, check.names = FALSE)
feature_cols <- colnames(data)[1:11]

# 2. Parameter Configuration
# ------------------------------------------------------------------------------
main_font <- "sans"
gap_size  <- 40      
h_scale   <- 1.8     
common_label_size <- 3.8 

# Geometric constraints for the central boxes
box_width_data      <- 0.8   
img_box_half_height <- 28    
desired_ratio       <- box_width_data / (img_box_half_height * 2)

# Color Palettes
grad_low  <- "#64B5F6" 
grad_high <- "#0D47A1"
col_func  <- colorRamp(c(grad_low, grad_high))

color_bot_tot  <- "#81D4FA"  
color_bot_npr  <- "#CE93D8"  
color_bot_prc  <- "#4DB6AC"  

bg_color_stroma <- "#F5F9FF" 
bg_color_tumor  <- "#FFF9FA" 
line_size       <- 0.25      

# 3. Data Reconstruction
# ------------------------------------------------------------------------------
plot_list <- list()

# Calculate reference means for normalization
all_ref_means <- sapply(feature_cols, function(f) mean(data[data$END1 == 2, f]) * 100)
min_ref <- min(all_ref_means)
max_ref <- max(all_ref_means)

get_grad_hex <- function(val) {
  norm_val <- (val - min_ref) / (max_ref - min_ref)
  rgb_vec  <- col_func(norm_val)
  rgb(rgb_vec[1], rgb_vec[2], rgb_vec[3], maxColorValue = 255)
}

top_bar_left     <- -0.35
top_bar_right    <- 0.35
gap_between_cols <- 0.03 

for (feat in feature_cols) {
  ref_vals   <- data[data$END1 == 2, feat]
  total_vals <- data[data$END1 == 1, feat]
  sub0_vals  <- data[data$END2 == 0, feat]
  sub1_vals  <- data[data$END2 == 1, feat]
  
  # Percentage calculations
  pct_ref   <- mean(ref_vals) * 100
  pct_total <- mean(total_vals) * 100
  pct_sub0  <- mean(sub0_vals) * 100
  pct_sub1  <- mean(sub1_vals) * 100
  
  # Statistical testing (Fisher's Exact Test)
  p_val <- fisher.test(matrix(c(sum(total_vals), length(total_vals)-sum(total_vals), 
                                sum(ref_vals), length(ref_vals)-sum(ref_vals)), 2))$p.value
  get_star <- function(p) ifelse(p<0.001, "***", ifelse(p<0.01, "**", ifelse(p<0.05, "*", "")))
  
  # Top segment: FH-dRCC
  plot_list[[length(plot_list)+1]] <- data.frame(
    Feature = feat, Group = "FH-dRCC", X_center = 0, 
    X_start = top_bar_left, X_end = top_bar_right, 
    Ymin = gap_size, Ymax = gap_size + pct_ref * h_scale,
    Label = sprintf("%.1f%%", pct_ref), Star = get_star(p_val), FillColor = get_grad_hex(pct_ref)
  )
  
  # Bottom segment: Total
  tot_end   <- -0.02
  plot_list[[length(plot_list)+1]] <- data.frame(
    Feature = feat, Group = "Total", 
    X_center = (top_bar_left + tot_end) / 2,
    X_start = top_bar_left, X_end = tot_end, 
    Ymin = -gap_size - pct_total * h_scale, Ymax = -gap_size,
    Label = sprintf("%.1f%%", pct_total), Star = "", FillColor = color_bot_tot
  )
  
  # Bottom segment: non-PRCC
  npr_start <- tot_end + gap_between_cols
  npr_end   <- 0.17
  plot_list[[length(plot_list)+1]] <- data.frame(
    Feature = feat, Group = "non-PRCC", 
    X_center = (npr_start + npr_end) / 2,
    X_start = npr_start, X_end = npr_end, 
    Ymin = -gap_size - pct_sub0 * h_scale, Ymax = -gap_size,
    Label = sprintf("%.1f%%", pct_sub0), Star = "", FillColor = color_bot_npr
  )
  
  # Bottom segment: PRCC
  prc_start <- npr_end + gap_between_cols
  plot_list[[length(plot_list)+1]] <- data.frame(
    Feature = feat, Group = "PRCC", 
    X_center = (prc_start + top_bar_right) / 2,
    X_start = prc_start, X_end = top_bar_right, 
    Ymin = -gap_size - pct_sub1 * h_scale, Ymax = -gap_size,
    Label = sprintf("%.1f%%", pct_sub1), Star = "", FillColor = color_bot_prc
  )
}

df <- do.call(rbind, plot_list)
df$Feature <- factor(df$Feature, levels = feature_cols)
img_box_data <- data.frame(Feature = feature_cols, x = 1:length(feature_cols))
y_max_limit <- max(df$Ymax) + 25 
y_min_limit <- min(df$Ymin) - 30 

# 4. Visualization
# ------------------------------------------------------------------------------
p <- ggplot() +
  # Background regions
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = bg_color_stroma) +
  annotate("rect", xmin = 4.5, xmax = 11.5, ymin = -Inf, ymax = Inf, fill = bg_color_tumor) +
  
  # Region headers
  annotate("text", x = 2.5, y = y_max_limit - 5, label = "Stroma characteristics", 
           size = 7, fontface = "bold", color = "black", family = main_font) + 
  annotate("text", x = 8, y = y_max_limit - 5, label = "Tumoral characteristics", 
           size = 7, fontface = "bold", color = "black", family = main_font) + 
  
  # Central image containers
  geom_rect(data = img_box_data,
            aes(xmin = x - box_width_data/2, xmax = x + box_width_data/2, 
                ymin = -img_box_half_height, ymax = img_box_half_height),
            fill = "white", color = "black", linewidth = 0.2) + 
  
  # Main bar charts
  geom_rect(data = df,
            aes(xmin = as.numeric(Feature) + X_start, xmax = as.numeric(Feature) + X_end,
                ymin = Ymin, ymax = Ymax, fill = FillColor),
            color = "black", linewidth = line_size) +
  
  # Labels and Significance stars
  geom_text(data = subset(df, Group == "FH-dRCC"),
            aes(x = as.numeric(Feature), y = Ymax + 2, label = Label, color = FillColor),
            vjust = -0.5, size = common_label_size, fontface = "bold") +
  geom_text(data = subset(df, Group == "FH-dRCC"),
            aes(x = as.numeric(Feature), y = Ymax + 8, label = Star),
            vjust = 0, size = 6, fontface = "bold", color = "black") +
  geom_text(data = df[df$Group != "FH-dRCC",],
            aes(x = as.numeric(Feature) + X_center, y = Ymin - 3, label = Label, color = FillColor),
            vjust = 1, size = common_label_size, fontface = "bold") +
  
  # Aspect Ratio and Coordinate locking
  coord_fixed(ratio = desired_ratio, xlim = c(0.5, 11.5), ylim = c(y_min_limit, y_max_limit), expand = FALSE) +
  
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(breaks = 1:length(feature_cols), labels = feature_cols) +
  
  # Axis Titles
  annotate("text", x = 0.2, y = gap_size + 30, label = "High scores (FH-dRCC) %", 
           angle = 90, size = 4.5, color = "black", fontface="bold") +
  annotate("text", x = 0.2, y = -gap_size - 30, label = "High scores (non-FH) %", 
           angle = 90, size = 4.5, color = "black", fontface="bold") +
  
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14, base_family = main_font) + 
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_blank(),
    plot.margin = margin(30, 10, 30, 10), 
    plot.background = element_rect(fill = "white", color = NA)
  )

# 5. Legend Construction
# ------------------------------------------------------------------------------
draw_key <- function(p, x, y, col, txt) {
  p + 
    annotate("rect", xmin = x, xmax = x + 0.5, ymin = y, ymax = y + 5, 
             fill = col, color = "black", linewidth = line_size) +
    annotate("text", x = x + 0.7, y = y + 2.5, label = txt, 
             hjust = 0, size = 3.5, color = "black")
}

legend_x <- 9.8
legend_base_y <- y_min_limit + 10

p <- p + annotate("text", x = legend_x, y = legend_base_y + 24, 
                  label = "non-FH-dRCC", hjust = 0, fontface = "bold", size = 4)

p <- draw_key(p, legend_x, legend_base_y + 16, color_bot_tot, "Total")
p <- draw_key(p, legend_x, legend_base_y + 8,  color_bot_npr, "non-PRCC")
p <- draw_key(p, legend_x, legend_base_y + 0,  color_bot_prc, "PRCC")

# 6. Final Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(file.path(results_dir, "Figure_4B_Analysis.pdf"), p, width = 16, height = 8)
message("Visualization saved to /results/ directory.")


# ==============================================================================
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
# Script Name: Figure_5_ROC_Analysis.R
# Description: Generate ROC curves for Figure 5.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here)

# Load data from the 'data' directory
data_path <- here("data", "Figure 5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
datasets <- unique(raw_data$Dataset)
raw_data$Dataset <- factor(raw_data$Dataset, levels = datasets)

# Define consistent color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 3. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list() 

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (AUC & CI) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  # Prepare annotation text
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Generate ROC Coordinates ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  # --- Step C: Create Individual Plot ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    # Draw curve with standardized linewidth
    geom_path(color = custom_colors[ds], linewidth = 0.8) +
    
    # Reference diagonal line
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 4, hjust = 0, lineheight = 1.2) +
    
    # Configure axes with reversed specificity and padding
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = ds) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[ds]] <- p
}

# 4. Final Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure5_ROC_Analysis.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Combined plot saved to: ", results_dir)


# ==============================================================================
# Script Name: Figure_5_Confusion_Matrices.R
# Description: Generate confusion matrices for datasets in Figure 5.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "caret", "gridExtra", "scales", "dplyr", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(caret)
library(gridExtra)
library(scales)
library(dplyr)
library(here)

# 2. Data Import
# ------------------------------------------------------------------------------
# Load data from the 'data' directory
data_path <- here("data", "Figure 5.csv")

if (!file.exists(data_path)) {
  stop("Data file not found! Please ensure 'Figure 5.csv' is in the /data/ directory.")
}

raw_data <- read.csv(data_path)

# 3. Data Preprocessing
# ------------------------------------------------------------------------------
# Ensure labels are factors with consistent levels
raw_data$true_labels <- factor(raw_data$true_labels, levels = c(0, 1))
raw_data$pred_labels <- factor(raw_data$pred_labels, levels = c(0, 1))

datasets <- unique(raw_data$Dataset)

# 4. Define Color Palettes
# ------------------------------------------------------------------------------
color_palettes <- list()
color_palettes[["Dataset2---Internal Testing"]] <- c("#FEFDFB", "#C1DEE6", "#79B5D0", "#1C93B5")
color_palettes[["Dataset3---Multicenter"]]       <- c("#FEFDFB", "#C5D5B8", "#97B786", "#3C7526")
color_palettes[["Dataset4---TCGA cohort"]]       <- c("#FEFDFB", "#FFE2BC", "#FFBE67", "#F2A312")

# 5. Plotting Function
# ------------------------------------------------------------------------------
plot_cm_square <- function(cm_table, dataset_name, palette_colors) {
  
  # Calculate percentages for heatmap intensity
  total_count <- sum(cm_table$Freq)
  cm_table$Freq_percentage <- cm_table$Freq / total_count
  
  # Generate Heatmap
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    # Tile layer with standardized linewidth
    geom_tile(color = "white", linewidth = 1.2) + 
    
    # Add count and percentage labels
    geom_text(aes(label = paste0(Freq, "\n", "(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 5) +
    
    scale_fill_gradientn(colors = palette_colors,
                         limits = c(0, 1),
                         labels = scales::percent,
                         name = NULL) +
    
    scale_x_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    scale_y_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    
    labs(title = dataset_name, x = "Predicted", y = "Actual") +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11, color = "black"),
      panel.grid = element_blank(), 
      legend.position = "right"
    ) +
    
    coord_fixed()
  
  return(p)
}

# 6. Generate Plots
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(raw_data, Dataset == ds)
  
  # Calculate Confusion Matrix using caret
  cm_obj <- confusionMatrix(sub_data$pred_labels, sub_data$true_labels)
  cm_df <- as.data.frame(cm_obj$table)
  
  # Assign palette (Fallback to grayscale if missing)
  current_palette <- color_palettes[[ds]]
  if (is.null(current_palette)) current_palette <- c("#FFFFFF", "#CCCCCC", "#666666", "#000000")
  
  plot_list[[ds]] <- plot_cm_square(cm_df, ds, current_palette)
}

# 7. Export Results
# ------------------------------------------------------------------------------
output_dir <- here("results")
if (!dir.exists(output_dir)) dir.create(output_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(output_dir, "Figure5_Confusion_Matrices.pdf"), 
       plot = combined_plot, width = 15, height = 5, device = "pdf")

message("Analysis complete. Output saved to /results/ directory.")


# ==============================================================================
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
# Script Name: Table_1_Performance_Analysis.R
# Description: Calculate comprehensive performance metrics (Sens, Spec, PPV, 
#              NPV, Brier Score, AUROC) including an "Overall" dataset.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)   # For ROC/AUC analysis
library(dplyr)  # For data manipulation
library(here)   # For relative path management

# Load data using relative path
data_path <- here("data", "Table 1.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preparation
# ------------------------------------------------------------------------------
# Create "Overall Dataset" by duplicating all entries
overall_data <- raw_data
overall_data$Dataset <- "Overall Dataset"

# Combine individual datasets with the overall summary
combined_data <- rbind(raw_data, overall_data)

# Set factor levels to control the order (Individual datasets followed by Overall)
orig_names <- unique(raw_data$Dataset)
combined_data$Dataset <- factor(combined_data$Dataset, levels = c(orig_names, "Overall Dataset"))

# 3. Metric Calculation Function
# ------------------------------------------------------------------------------
calculate_metrics_final <- function(group_data) {
  
  # Ensure labels and predictions are numeric
  actual <- as.numeric(as.character(group_data$true_labels))
  predicted <- as.numeric(as.character(group_data$pred_labels)) 
  
  # --- Part 1: Sensitivity, Specificity, PPV, NPV ---
  # Construct confusion matrix
  cm <- table(factor(actual, levels = c(0, 1)), 
              factor(predicted, levels = c(0, 1)))
  
  TP <- cm["1", "1"]; FN <- cm["1", "0"]
  TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  # Helper for calculating proportions with 95% CI
  get_prop_ci <- function(k, n) {
    if (n > 0) {
      res <- prop.test(k, n, conf.level = 0.95)
      return(c(val = k/n, lower = res$conf.int[1], upper = res$conf.int[2]))
    } else {
      return(c(val = 0, lower = 0, upper = 0))
    }
  }
  
  sens <- get_prop_ci(TP, TP + FN)
  spec <- get_prop_ci(TN, TN + FP)
  ppv  <- get_prop_ci(TP, TP + FP)
  npv  <- get_prop_ci(TN, TN + FN)
  
  # --- Part 2: AUROC (DeLong Method) ---
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- ci.auc(roc_obj, conf.level = 0.95, method = "delong")
  
  # --- Part 3: Brier Score (Bootstrap CI) ---
  brier_val <- mean((predicted - actual)^2)
  
  set.seed(123) # Ensure reproducibility for bootstrap
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), length(actual), replace = TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs = c(0.025, 0.975))
  
  # --- Part 4: Formatting Result ---
  fmt <- function(val, low, high) {
    sprintf("%.3f (%.3f-%.3f)", val, low, high)
  }
  
  result <- data.frame(
    Dataset = unique(group_data$Dataset),
    Sensitivity = fmt(sens['val'], sens['lower'], sens['upper']),
    Specificity = fmt(spec['val'], spec['lower'], spec['upper']),
    PPV         = fmt(ppv['val'],  ppv['lower'],  ppv['upper']),
    NPV         = fmt(npv['val'],  npv['lower'],  npv['upper']),
    Brier_Score = fmt(brier_val, brier_ci[1], brier_ci[2]),
    AUROC       = fmt(auc_val, auc_ci[1], auc_ci[3]),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# 4. Execution and Export
# ------------------------------------------------------------------------------
final_results <- combined_data %>%
  group_by(Dataset) %>%
  do(calculate_metrics_final(.)) %>%
  ungroup()

# Export results to the 'results' directory
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_results, file.path(results_dir, "Table_1_Performance_Metrics.csv"), row.names = FALSE)

# Display results in console
print(final_results)
message("Calculation complete. Results exported to /results/ directory.")


# ==============================================================================
# Script Name: Table_S3_Patient_Analysis.R
# Description: Calculate aggregated patient-level performance metrics.
#              Rule: Max pred_label per patient is used for diagnosis.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S3.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure labels are numeric
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# 2. Patient-Level Aggregation
# ------------------------------------------------------------------------------
# Robust identification of Patient ID column (matching "Pts")
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (is.na(pts_col)) stop("Patient ID column (Pts) not found.")
colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Aggregation: Assign patient-level diagnosis based on the maximum label per patient
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>% 
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Add Overall Summary
# ------------------------------------------------------------------------------
overall_data <- patient_data
overall_data$Dataset <- "Overall Dataset"

combined_data <- rbind(patient_data, overall_data)

# Control dataset order: Individual sets first, followed by Overall summary
orig_names <- unique(patient_data$Dataset)
combined_data$Dataset <- factor(combined_data$Dataset, levels = c(orig_names, "Overall Dataset"))

# 4. Metrics Calculation Function
# ------------------------------------------------------------------------------
calc_metrics <- function(df) {
  
  actual <- df$true_labels
  predicted <- df$pred_labels
  
  # --- Contingency Table Metrics (Sens, Spec, PPV, NPV) ---
  cm <- table(factor(actual, levels=c(0,1)), factor(predicted, levels=c(0,1)))
  TP <- cm["1","1"]; FN <- cm["1","0"]; TN <- cm["0","0"]; FP <- cm["0","1"]
  
  # Proportion test for confidence intervals
  get_ci <- function(k, n) {
    if(n > 0) {
      res <- prop.test(k, n, conf.level=0.95)
      c(res$estimate, res$conf.int)
    } else c(0, 0, 0)
  }
  
  sens <- get_ci(TP, TP+FN)
  spec <- get_ci(TN, TN+FP)
  ppv  <- get_ci(TP, TP+FP)
  npv  <- get_ci(TN, TN+FN)
  
  # --- AUROC (DeLong Method) ---
  roc_obj <- roc(actual, predicted, quiet=TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  # Using tryCatch to handle edge cases where DeLong might fail
  auc_ci <- tryCatch(ci.auc(roc_obj, method="delong"), error=function(e) c(NA,NA,NA))
  
  # --- Brier Score (Bootstrap) ---
  brier_val <- mean((predicted - actual)^2)
  set.seed(123) # Reproducibility
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), replace=TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs=c(0.025, 0.975))
  
  # --- Formatting Results ---
  fmt <- function(v) sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
  
  data.frame(
    Dataset = unique(df$Dataset),
    Sensitivity = fmt(sens),
    Specificity = fmt(spec),
    PPV = fmt(ppv),
    NPV = fmt(npv),
    Brier_Score = sprintf("%.3f (%.3f-%.3f)", brier_val, brier_ci[1], brier_ci[2]),
    AUROC = if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA",
    stringsAsFactors = FALSE
  )
}

# 5. Execution and Export
# ------------------------------------------------------------------------------
final_results <- combined_data %>%
  group_by(Dataset) %>%
  do(calc_metrics(.)) %>%
  ungroup()

# Create results directory if needed
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_results, file.path(results_dir, "Table_S3_Patient_Metrics.csv"), row.names = FALSE)

# Display Results
print(final_results)
message("Patient-level metrics calculation complete. Results saved to /results/.")


# ==============================================================================
# Script Name: Table_S5_Subgroup_Analysis.R
# Description: Subgroup performance metrics for Table S5.
#              Calculates Sens, Spec, PPV, NPV, Brier Score, and AUROC
#              using binary 'pred_labels'.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure numeric format for labels
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# 2. Metric Calculation Function
# ------------------------------------------------------------------------------
calculate_metrics <- function(df) {
  
  actual    <- df$true_labels
  predicted <- df$pred_labels 
  
  # --- Contingency Table Metrics (Sens, Spec, PPV, NPV) ---
  cm <- table(factor(actual, levels = c(0, 1)), factor(predicted, levels = c(0, 1)))
  TP <- cm["1", "1"]; FN <- cm["1", "0"]; TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  # Proportion test for 95% Confidence Intervals
  get_ci <- function(k, n) {
    if (n > 0) {
      res <- prop.test(k, n, conf.level = 0.95)
      c(res$estimate, res$conf.int)
    } else c(NA, NA, NA)
  }
  
  sens <- get_ci(TP, TP + FN)
  spec <- get_ci(TN, TN + FP)
  ppv  <- get_ci(TP, TP + FP)
  npv  <- get_ci(TN, TN + FN)
  
  # --- AUROC (DeLong Method using Binary Labels) ---
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) c(NA, NA, NA))
  
  # --- Brier Score (Bootstrap) ---
  brier_val <- mean((predicted - actual)^2)
  
  set.seed(123) # For reproducibility
  brier_boot <- replicate(1000, {
    idx <- sample(seq_along(actual), replace = TRUE)
    mean((predicted[idx] - actual[idx])^2)
  })
  brier_ci <- quantile(brier_boot, probs = c(0.025, 0.975))
  
  # --- Formatting Results ---
  fmt <- function(v) {
    if(all(is.na(v))) return("NA")
    sprintf("%.3f (%.3f-%.3f)", v[1], v[2], v[3])
  }
  
  data.frame(
    Group       = unique(df$GROUP),
    Sensitivity = fmt(sens),
    Specificity = fmt(spec),
    PPV         = fmt(ppv),
    NPV         = fmt(npv),
    Brier_Score = sprintf("%.3f (%.3f-%.3f)", brier_val, brier_ci[1], brier_ci[2]),
    AUROC       = if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA",
    stringsAsFactors = FALSE
  )
}

# 3. Execution and Export
# ------------------------------------------------------------------------------
final_results <- raw_data %>%
  group_by(GROUP) %>%
  do(calculate_metrics(.)) %>%
  ungroup()

# Create results directory if needed
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Save to CSV
write.csv(final_results, file.path(results_dir, "Table_S5_Metrics_Binary.csv"), row.names = FALSE)

# Display Results
print(final_results)
message("Subgroup analysis complete. Results saved to /results/.")


# ==============================================================================
# Script Name: Table_S6_Reader_Comparison.R
# Description: 
#   1. Compares multiple predictors against a Reference.
#   2. Performs McNemar's test (Sens/Spec), Bootstrap (PPV/NPV), and DeLong (AUC).
#   3. Generates summary rows for Junior, Senior, and All pathologists.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
required_packages <- c("pROC", "dplyr", "stringr", "here")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(pROC)
library(dplyr)
library(stringr)
library(here)

# Load data using relative path
data_path <- here("data", "Table S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure all columns are treated as numeric for analysis
raw_data <- raw_data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

# 2. Helper Functions
# ------------------------------------------------------------------------------
# Format proportion with 95% Confidence Interval
fmt_ci <- function(k, n) {
  if (n > 0) {
    res <- prop.test(k, n, conf.level = 0.95)
    sprintf("%.3f (%.3f-%.3f)", k/n, res$conf.int[1], res$conf.int[2])
  } else "NA"
}

# Format P-values for publication
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# 3. Core Analysis Engine
# ------------------------------------------------------------------------------
analyze_comparator <- function(target, ref_pred, comp_pred, col_name, is_ref = FALSE) {
  
  # --- Metric Calculation ---
  cm <- table(factor(target, levels = c(0, 1)), factor(comp_pred, levels = c(0, 1)))
  TP <- cm["1", "1"]; FN <- cm["1", "0"]
  TN <- cm["0", "0"]; FP <- cm["0", "1"]
  
  val_sens <- fmt_ci(TP, TP + FN)
  val_spec <- fmt_ci(TN, TN + FP)
  val_ppv  <- fmt_ci(TP, TP + FP)
  val_npv  <- fmt_ci(TN, TN + FN)
  
  roc_comp <- roc(target, comp_pred, quiet = TRUE)
  auc_val  <- as.numeric(auc(roc_comp))
  auc_ci   <- tryCatch(ci.auc(roc_comp, method = "delong"), error = function(e) c(NA, NA, NA))
  val_auc  <- if(!is.na(auc_ci[1])) sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci[1], auc_ci[3]) else "NA"
  
  # --- Statistical Comparison against Reference ---
  p_sens <- "Ref."; p_spec <- "Ref."; p_ppv <- "Ref."; p_npv <- "Ref."; p_auc <- "Ref."
  
  if (!is_ref) {
    # McNemar's Test for Sensitivity and Specificity
    idx_pos <- which(target == 1)
    if(length(idx_pos) > 0) {
      tbl_sens <- table(factor(ref_pred[idx_pos], levels=c(0,1)), factor(comp_pred[idx_pos], levels=c(0,1)))
      p_sens <- fmt_p(tryCatch(mcnemar.test(tbl_sens)$p.value, error=function(e) NA))
    }
    idx_neg <- which(target == 0)
    if(length(idx_neg) > 0) {
      tbl_spec <- table(factor(ref_pred[idx_neg], levels=c(0,1)), factor(comp_pred[idx_neg], levels=c(0,1)))
      p_spec <- fmt_p(tryCatch(mcnemar.test(tbl_spec)$p.value, error=function(e) NA))
    }
    
    # DeLong's Test for AUC Comparison
    roc_ref <- roc(target, ref_pred, quiet = TRUE)
    test_res <- tryCatch(roc.test(roc_ref, roc_comp, method = "delong"), error=function(e) NULL)
    if(!is.null(test_res)) p_auc <- fmt_p(test_res$p.value)
    
    # Bootstrap for PPV and NPV comparison (1000 iterations)
    calc_diff <- function(data, idx) {
      d <- data[idx, ]
      cm_r <- table(factor(d$t, 0:1), factor(d$r, 0:1))
      ppv_r <- cm_r[2,2]/sum(cm_r[,2]); npv_r <- cm_r[1,1]/sum(cm_r[,1])
      cm_c <- table(factor(d$t, 0:1), factor(d$c, 0:1))
      ppv_c <- cm_c[2,2]/sum(cm_c[,2]); npv_c <- cm_c[1,1]/sum(cm_c[,1])
      c(ppv_c - ppv_r, npv_c - npv_r)
    }
    
    set.seed(123)
    boot_dat <- data.frame(t=target, r=ref_pred, c=comp_pred)
    boot_diffs <- replicate(1000, calc_diff(boot_dat, sample(nrow(boot_dat), replace=TRUE)))
    
    get_p <- function(vec) {
      p <- 2 * min(mean(vec > 0, na.rm=TRUE), mean(vec < 0, na.rm=TRUE))
      fmt_p(min(p, 1.0))
    }
    p_ppv <- get_p(boot_diffs[1,]); p_npv <- get_p(boot_diffs[2,])
  }
  
  data.frame(Reader = col_name, 
             Sensitivity = val_sens, P_Sens = p_sens,
             Specificity = val_spec, P_Spec = p_spec,
             PPV = val_ppv, P_PPV = p_ppv,
             NPV = val_npv, P_NPV = p_npv,
             AUROC = val_auc, P_AUC = p_auc, stringsAsFactors = FALSE)
}

# 4. Main Table Generation
# ------------------------------------------------------------------------------
target_col <- raw_data[[1]]
ref_col    <- raw_data[[2]]
col_names  <- colnames(raw_data)

message("Processing individual columns...")
results_list <- lapply(2:ncol(raw_data), function(i) {
  analyze_comparator(target_col, ref_col, raw_data[[i]], col_names[i], is_ref = (i == 2))
})

final_table <- do.call(rbind, results_list)

# 5. Summary Rows (Averages)
# ------------------------------------------------------------------------------
# Extract numeric value from formatted string "0.xxx (0.xxx-0.xxx)"
get_val <- function(str_val) {
  if (str_val == "NA" || is.na(str_val)) return(NA)
  as.numeric(str_split(str_val, " ")[[1]][1])
}

# Calculate average row based on specified indices
calc_avg_row <- function(df, row_indices, label_name) {
  metric_cols <- c("Sensitivity", "Specificity", "PPV", "NPV", "AUROC")
  new_row <- list(Reader = label_name)
  
  for (col in colnames(df)) {
    if (col == "Reader") next
    if (col %in% metric_cols) {
      vals <- sapply(df[row_indices, col], get_val)
      new_row[[col]] <- sprintf("%.3f", mean(vals, na.rm = TRUE))
    } else {
      new_row[[col]] <- "-" # No statistical test for average rows
    }
  }
  return(as.data.frame(new_row, stringsAsFactors = FALSE))
}

# Define grouping indices (Column 3-6: Junior, Column 7-9: Senior)
idx_junior <- (3:6) - 1
idx_senior <- (7:9) - 1
idx_all    <- (3:9) - 1

avg_junior <- calc_avg_row(final_table, idx_junior, "Average (Junior pathologists)")
avg_senior <- calc_avg_row(final_table, idx_senior, "Average (Senior pathologists)")
avg_all    <- calc_avg_row(final_table, idx_all,    "Average (All pathologists)")

final_table_with_avg <- rbind(final_table, avg_junior, avg_senior, avg_all)

# 6. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_table_with_avg, file.path(results_dir, "Table_S6_Comparison_Results.csv"), row.names = FALSE)

print(final_table_with_avg)
message("Analysis complete. Results exported to /results/.")


# ==============================================================================
# Script Name: Table_S7_Paired_Analysis.R
# Description: Paired performance analysis by pathologist group.
#              Compares "Modified strategy" vs "Without AI".
#              Metrics: Sens, Spec, PPV, NPV, AUROC, Brier Score.
#              Tests: McNemar, DeLong, and Bootstrap for paired differences.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(pROC)
library(dplyr)
library(stringr)
library(here)

# Load data using relative path
# check.names = FALSE ensures special characters in headers are preserved
data_path <- here("data", "Table S7.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path, check.names = FALSE)

# Ensure numeric format across all diagnostic columns
raw_data <- raw_data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

# 2. Helper Functions for Formatting and Stats
# ------------------------------------------------------------------------------
fmt_ci <- function(k, n) {
  if (n > 0) {
    res <- prop.test(k, n, conf.level = 0.95)
    sprintf("%.3f (%.3f-%.3f)", k/n, res$conf.int[1], res$conf.int[2])
  } else "NA"
}

fmt_val_ci <- function(val, lower, upper) {
  sprintf("%.3f (%.3f-%.3f)", val, lower, upper)
}

fmt_p <- function(p) {
  if (is.na(p) || is.nan(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# 3. Paired Group Analysis Function
# ------------------------------------------------------------------------------
analyze_paired_group <- function(target, ref_pred, comp_pred, group_name) {
  
  # --- A. Metrics for "Without AI" (Reference) ---
  cm_r <- table(factor(target, levels=0:1), factor(ref_pred, levels=0:1))
  sens_r <- fmt_ci(cm_r[2,2], sum(cm_r[2,])); spec_r <- fmt_ci(cm_r[1,1], sum(cm_r[1,]))
  ppv_r  <- fmt_ci(cm_r[2,2], sum(cm_r[,2])); npv_r  <- fmt_ci(cm_r[1,1], sum(cm_r[,1]))
  
  roc_r <- roc(target, ref_pred, quiet=TRUE)
  auc_r_val <- as.numeric(auc(roc_r))
  auc_r_ci  <- tryCatch(ci.auc(roc_r, method="delong"), error=function(e) c(NA,NA,NA))
  auc_r <- fmt_val_ci(auc_r_val, auc_r_ci[1], auc_r_ci[3])
  
  brier_r_val <- mean((ref_pred - target)^2)
  set.seed(123)
  b_boot_r <- replicate(1000, { idx<-sample(length(target), replace=T); mean((ref_pred[idx]-target[idx])^2) })
  brier_r_ci <- quantile(b_boot_r, c(0.025, 0.975))
  brier_r <- fmt_val_ci(brier_r_val, brier_r_ci[1], brier_r_ci[2])
  
  # --- B. Metrics for "Modified strategy" (Comparator) ---
  cm_c <- table(factor(target, levels=0:1), factor(comp_pred, levels=0:1))
  sens_c <- fmt_ci(cm_c[2,2], sum(cm_c[2,])); spec_c <- fmt_ci(cm_c[1,1], sum(cm_c[1,]))
  ppv_c  <- fmt_ci(cm_c[2,2], sum(cm_c[,2])); npv_c  <- fmt_ci(cm_c[1,1], sum(cm_c[,1]))
  
  roc_c <- roc(target, comp_pred, quiet=TRUE)
  auc_c_val <- as.numeric(auc(roc_c))
  auc_c_ci  <- tryCatch(ci.auc(roc_c, method="delong"), error=function(e) c(NA,NA,NA))
  auc_c <- fmt_val_ci(auc_c_val, auc_c_ci[1], auc_c_ci[3])
  
  brier_c_val <- mean((comp_pred - target)^2)
  set.seed(123)
  b_boot_c <- replicate(1000, { idx<-sample(length(target), replace=T); mean((comp_pred[idx]-target[idx])^2) })
  brier_c_ci <- quantile(b_boot_c, c(0.025, 0.975))
  brier_c <- fmt_val_ci(brier_c_val, brier_c_ci[1], brier_c_ci[2])
  
  # --- C. Statistical Comparisons ---
  idx_pos <- which(target==1); idx_neg <- which(target==0)
  
  p_sens <- tryCatch(mcnemar.test(table(factor(ref_pred[idx_pos],0:1), factor(comp_pred[idx_pos],0:1)))$p.value, error=function(e) NA)
  p_spec <- tryCatch(mcnemar.test(table(factor(ref_pred[idx_neg],0:1), factor(comp_pred[idx_neg],0:1)))$p.value, error=function(e) NA)
  p_auc <- tryCatch(roc.test(roc_r, roc_c, method="delong")$p.value, error=function(e) NA)
  
  # Bootstrap for PPV, NPV, and Brier Score differences
  calc_diffs <- function(d, idx) {
    curr <- d[idx, ]
    tb_r <- table(factor(curr$t,0:1), factor(curr$r,0:1))
    ppv_r <- tb_r[2,2]/sum(tb_r[,2]); npv_r <- tb_r[1,1]/sum(tb_r[,1]); br_r <- mean((curr$r - curr$t)^2)
    tb_c <- table(factor(curr$t,0:1), factor(curr$c,0:1))
    ppv_c <- tb_c[2,2]/sum(tb_c[,2]); npv_c <- tb_c[1,1]/sum(tb_c[,1]); br_c <- mean((curr$c - curr$t)^2)
    c(ppv_c - ppv_r, npv_c - npv_r, br_c - br_r)
  }
  
  boot_df <- data.frame(t=target, r=ref_pred, c=comp_pred)
  set.seed(123)
  boot_res <- replicate(1000, calc_diffs(boot_df, sample(nrow(boot_df), replace=TRUE)))
  
  get_p <- function(vec) {
    p <- 2 * min(mean(vec > 0, na.rm=TRUE), mean(vec < 0, na.rm=TRUE))
    fmt_p(min(p, 1.0))
  }
  p_ppv <- get_p(boot_res[1,]); p_npv <- get_p(boot_res[2,]); p_brier <- get_p(boot_res[3,])
  
  # --- D. Formatting Output Table ---
  row_comp <- data.frame(
    Group = group_name, Strategy = "Modified strategy",
    Sensitivity = sens_c, P_Sens = fmt_p(p_sens),
    Specificity = spec_c, P_Spec = fmt_p(p_spec),
    PPV = ppv_c, P_PPV = p_ppv,
    NPV = npv_c, P_NPV = p_npv,
    AUROC = auc_c, P_AUC = fmt_p(p_auc),
    Brier_Score = brier_c, P_Brier = p_brier,
    stringsAsFactors = FALSE
  )
  
  row_ref <- data.frame(
    Group = group_name, Strategy = "Without AI",
    Sensitivity = sens_r, P_Sens = "Ref.",
    Specificity = spec_r, P_Spec = "Ref.",
    PPV = ppv_r, P_PPV = "Ref.",
    NPV = npv_r, P_NPV = "Ref.",
    AUROC = auc_r, P_AUC = "Ref.",
    Brier_Score = brier_r, P_Brier = "Ref.",
    stringsAsFactors = FALSE
  )
  
  return(rbind(row_comp, row_ref))
}

# 4. Processing Loop
# ------------------------------------------------------------------------------
target_col <- raw_data[[1]]
all_colnames <- colnames(raw_data)[-1]
groups <- unique(sapply(str_split(all_colnames, "_"), `[`, 1))

results_list <- list()

message("Executing paired analysis...")

for (grp in groups) {
  cols_in_grp <- all_colnames[startsWith(all_colnames, paste0(grp, "_"))]
  col_ref <- cols_in_grp[grep("Without AI", cols_in_grp)]
  col_comp <- cols_in_grp[grep("Modified strategy", cols_in_grp)]
  
  if (length(col_ref) >= 1 && length(col_comp) >= 1) {
    results_list[[grp]] <- analyze_paired_group(target_col, raw_data[[col_ref[1]]], raw_data[[col_comp[1]]], grp)
  }
}

final_output <- do.call(rbind, results_list)

# 5. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

write.csv(final_output, file.path(results_dir, "Table_S7_Paired_Comparison.csv"), row.names = FALSE)

print(final_output)
message("Analysis complete. Results exported to /results/.")


# ==============================================================================
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
# Script Name: Figure_S4_Threshold_Analysis.R
# Description: Multi-dataset confusion matrix analysis across varying cutoffs.
#              - Rows: Datasets; Columns: Probability Cutoffs (0.3 - 0.8).
#              - Independent color scales and legends per dataset.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(caret)
library(gridExtra)
library(grid)
library(scales) 
library(cowplot)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S4.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

data <- read.csv(data_path)
data$true_labels <- as.numeric(as.character(data$true_labels))
data$pred_probs  <- as.numeric(as.character(data$pred_probs))

datasets <- unique(data$Dataset)
cutoffs <- seq(0.3, 0.8, by = 0.1)

# 2. Color Palette Configuration
# ------------------------------------------------------------------------------
color_palettes <- list(
  "Dataset1---Training"         = c("#FEFDFB","#D4D1E6","#A29CCC","#7673A6"),
  "Dataset2---Internal Testing" = c("#FEFDFB","#C1DEE6","#79B5D0","#1C93B5"),
  "Dataset3---Multicenter"      = c("#FEFDFB","#C5D5B8","#97B786","#3C7526"),
  "Dataset4---TCGA cohort"      = c("#FEFDFB","#FFE2BC","#FFBE67","#F2A312")
)

# 3. Plotting Functions
# ------------------------------------------------------------------------------
# Function to generate individual confusion matrix tile
plot_confusion_matrix <- function(cm, cutoff, palette_colors, is_first_col, is_first_row) {
  cm_table <- as.data.frame(cm$table)
  cm_table$Freq_percentage <- cm_table$Freq / sum(cm_table$Freq)
  
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    geom_tile(color = "white", linewidth = 0.5) + 
    geom_text(aes(label = paste0(Freq, "\n(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 3) + 
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), name = NULL) +
    scale_x_discrete(labels = c("0" = "non-FH-dRCC", "1" = "FH-dRCC")) + 
    scale_y_discrete(labels = c("0" = "non-FH-dRCC", "1" = "FH-dRCC")) + 
    labs(
      title = if(is_first_row) paste0("Cutoff: ", cutoff) else NULL, 
      x = "Predicted", 
      y = if(is_first_col) "Actual" else NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = if(is_first_col) element_text(size = 9, face = "bold") else element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ) +
    coord_fixed()
  
  return(p)
}

# Function to extract colorbar legend for each row
get_row_legend <- function(palette_colors) {
  p <- ggplot(data.frame(x=1, y=1, val=0.5), aes(x, y, fill=val)) +
    geom_tile() +
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), 
                         labels = scales::percent, name = "Frequency") +
    guides(fill = guide_colorbar(barheight = unit(3.5, "cm"), barwidth = unit(0.4, "cm"),
                                 title.position = "top", title.hjust = 0.5)) +
    theme_minimal() +
    theme(legend.position = "right", 
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7))
  return(get_legend(p))
}

# 4. Processing Loop
# ------------------------------------------------------------------------------
plot_list <- list()
legend_list <- list() 
row_idx <- 0

for (ds in datasets) {
  row_idx <- row_idx + 1
  sub_data <- subset(data, Dataset == ds)
  current_palette <- color_palettes[[ds]] %||% c("#FFFFFF", "#CCCCCC", "#000000")
  
  legend_list[[ds]] <- get_row_legend(current_palette)
  
  col_idx <- 0
  for (cutoff in cutoffs) {
    col_idx <- col_idx + 1
    pred_factor <- factor(ifelse(sub_data$pred_probs >= cutoff, 1, 0), levels = c(0, 1))
    true_factor <- factor(sub_data$true_labels, levels = c(0, 1))
    
    cm <- confusionMatrix(pred_factor, true_factor)
    p <- plot_confusion_matrix(cm, cutoff, current_palette, (col_idx == 1), (row_idx == 1))
    plot_list[[paste0(ds, "_", cutoff)]] <- p
  }
}

# 5. Assembly and Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Construct sub-grids
main_grid <- arrangeGrob(grobs = plot_list, nrow = length(datasets), ncol = length(cutoffs))
legend_col <- arrangeGrob(grobs = legend_list, ncol = 1)
row_labels <- arrangeGrob(grobs = lapply(datasets, function(name) {
  textGrob(strsplit(as.character(name), "---")[[1]][1], rot = -90, gp = gpar(fontsize = 12, fontface = "bold"))
}), ncol = 1)

# Final Arrangement
final_plot <- grid.arrange(main_grid, legend_col, row_labels, ncol = 3, widths = c(10, 2, 0.5))

ggsave(filename = file.path(results_dir, "Figure_S4_Threshold_Analysis.pdf"), 
       plot = final_plot, width = length(cutoffs) * 2.2 + 2.5, height = length(datasets) * 2.4, device = "pdf")

message("Analysis complete. Results saved to /results/.")


# ==============================================================================
# Script Name: Figure_S5_DCA_and_PR_Analysis.R
# Description: Generates Decision Curves and Precision-Recall Curves 
#              across multiple cohorts to evaluate clinical utility.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggthemes)
library(rmda)
library(pROC)
library(PRROC)
library(here)

# Define dataset and output paths
data_path <- here("data", "Figure S5.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# 2. Configuration
# ------------------------------------------------------------------------------
group_col <- "Dataset" 
my_colors <- c("#7AB5D1", "#678655", "#F0A72C")

# Load Data
demo_data <- read.table(data_path, sep=",", header = TRUE)
if (!group_col %in% names(demo_data)) {
  stop(paste0("Error: Column '", group_col, "' not found."))
}

cohorts <- unique(demo_data[[group_col]])

# 3. PDF Export Initialization
# ------------------------------------------------------------------------------
output_filename <- file.path(results_dir, "Figure_S5_DCA_PR_Analysis.pdf")
message("Generating combined analysis PDF...")

pdf(output_filename, width = 10, height = 15) 

# Layout: 3 Rows, 2 Columns with outer margins
par(mfrow = c(3, 2), oma = c(2, 4, 1, 1))

# 4. Processing Loop
# ------------------------------------------------------------------------------
for (i in seq_along(cohorts)) {
  
  cohort_name <- cohorts[i]
  current_color <- my_colors[(i - 1) %% length(my_colors) + 1]
  sub_data <- demo_data[demo_data[[group_col]] == cohort_name, ]
  
  # --- Part A: Decision Curve Analysis (DCA) ---
  # Adjust margins for axis labels
  par(mar = c(7, 5, 3, 2)) 
  
  dca_model <- decision_curve(formula = true_labels ~ pred_probs, 
                              data = sub_data, family = binomial, 
                              thresholds = seq(0, 1, by = 0.01), policy = "opt-in")
  dca_model <- na.omit(dca_model)
  
  plot_decision_curve(dca_model,
                      curve.names = "Prediction Model",
                      col = current_color,
                      lwd = 3, lty = 1,
                      xlab = "Threshold Probability",
                      ylab = "Net Benefit",
                      cex.lab = 1.6, cex.axis = 1.4, cex.main = 1.8,
                      legend.position = "none",
                      cost.benefit.axis = TRUE)
  
  title("Decision Curve", col.main = current_color, font.main = 4, line = 1, cex.main = 1.8)
  
  # Add Dataset Label to the left margin
  mtext(cohort_name, side = 2, line = 5, cex = 1.5, font = 2, col = "black", las = 0)
  
  # --- Part B: Precision-Recall (PR) Curve ---
  par(mar = c(7, 5, 3, 2)) 
  
  true_labels <- sub_data$true_labels
  pred_probs <- sub_data$pred_probs
  
  pr <- pr.curve(scores.class0 = pred_probs[true_labels == 1],
                 scores.class1 = pred_probs[true_labels == 0], curve = TRUE)
  
  # Calculate ROC statistics for annotation
  roc_obj <- roc(true_labels, pred_probs, quiet = TRUE)
  auc_val <- auc(roc_obj)
  ci_val <- ci.auc(roc_obj)
  
  plot(pr,
       xlab = "Recall",
       ylab = "Precision",
       col = current_color,
       lwd = 3,
       cex.lab = 1.6, cex.axis = 1.4, cex.main = 1.8)
  
  # Annotate with AUC and Confidence Intervals
  text(x = 0, y = 0.05, 
       labels = sprintf("AUC: %.3f (95%% CI: %.3f-%.3f)", auc_val, ci_val[1], ci_val[3]),
       col = "black", cex = 1.3, adj = 0)
  
  box(lwd = 2)
  title("Recall Curve", col.main = current_color, font.main = 4, line = 1, cex.main = 1.8)
}

dev.off()
message("Analysis complete. PDF saved to: ", output_filename)


# ==============================================================================
# Script Name: Figure_S6_Patient_ROC.R
# Description: Generate Patient-Level ROC curves for Figure S6.
#              - Aggregation: Maximum probability and label per Patient ID.
#              - Curve Shape: Derived from 'pred_probs'.
#              - AUC Statistics: Derived from 'pred_labels'.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing
# ------------------------------------------------------------------------------
# Robustly identify Patient ID column
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (!is.na(pts_col)) colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Ensure numeric types
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_probs  <- as.numeric(as.character(raw_data$pred_probs))

# Identify prediction labels column (defaults to 5th column if name not found)
if(!"pred_labels" %in% colnames(raw_data)){
  colnames(raw_data)[5] <- "pred_labels"
}
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# Patient-Level Aggregation: Taking the maximum values
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>%
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_probs  = max(pred_probs, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

datasets <- unique(patient_data$Dataset)
patient_data$Dataset <- factor(patient_data$Dataset, levels = datasets)

# Define color palette
custom_colors <- c("#74B1CF", "#5F804D", "#F0A221")
names(custom_colors) <- datasets

# 3. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(patient_data, Dataset == ds)
  
  # --- Step A: Statistical Computation (Based on binary labels) ---
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  ci_result <- ci.auc(roc_stat_obj, method = "delong")
  
  label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                      "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  
  # --- Step B: Coordinate Extraction (Based on probabilities for smooth curve) ---
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Dataset = ds
  )
  
  current_color <- custom_colors[as.character(ds)]
  if(is.null(current_color)) current_color <- "black"
  
  # --- Step C: Create Individual Plot ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = current_color, linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 4, hjust = 0, lineheight = 1.2) +
    
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = ds) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[as.character(ds)]] <- p
}

# 4. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure_S6_Patient_ROC.pdf"), 
       plot = combined_plot, width = 13, height = 5, device = "pdf")

message("Analysis complete. Results exported to /results/.")


# ==============================================================================
# Script Name: Figure_S6_Patient_Confusion_Matrices.R
# Description: Generate Patient-Level Square Confusion Matrices for Figure S6.
#              Aggregation Rule: Max pred_labels per Patient ID.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(caret)
library(gridExtra)
library(scales)
library(dplyr)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S6.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# 2. Data Preprocessing (Aggregation)
# ------------------------------------------------------------------------------
# Identify Patient ID column (matching "Pts")
pts_col <- colnames(raw_data)[grep("Pts", colnames(raw_data), ignore.case = TRUE)][1]
if (!is.na(pts_col)) colnames(raw_data)[which(colnames(raw_data) == pts_col)] <- "Pts.ID"

# Ensure numeric types
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
if(!"pred_labels" %in% colnames(raw_data)) stop("Column 'pred_labels' not found.")
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))

# Aggregate: Take MAX pred_label for each patient (Patient-level diagnosis)
patient_data <- raw_data %>%
  group_by(Pts.ID, Dataset) %>%
  summarise(
    true_labels = max(true_labels, na.rm = TRUE),
    pred_labels = max(pred_labels, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to factors for Confusion Matrix
patient_data$true_labels <- factor(patient_data$true_labels, levels = c(0, 1))
patient_data$pred_labels <- factor(patient_data$pred_labels, levels = c(0, 1))

datasets <- unique(patient_data$Dataset)

# 3. Define Color Palettes
# ------------------------------------------------------------------------------
color_palettes <- list(
  "Dataset2---Internal Testing" = c("#FEFDFB", "#C1DEE6", "#79B5D0", "#1C93B5"),
  "Dataset3---Multicenter"      = c("#FEFDFB", "#C5D5B8", "#97B786", "#3C7526"),
  "Dataset4---TCGA cohort"      = c("#FEFDFB", "#FFE2BC", "#FFBE67", "#F2A312")
)

# 4. Plotting Function
# ------------------------------------------------------------------------------
plot_cm_square <- function(cm_table, dataset_name, palette_colors) {
  
  # Calculate percentages
  total_count <- sum(cm_table$Freq)
  cm_table$Freq_percentage <- cm_table$Freq / total_count
  
  # Generate Heatmap
  p <- ggplot(cm_table, aes(x = Prediction, y = Reference, fill = Freq_percentage)) +
    geom_tile(color = "white", linewidth = 1.2) + 
    geom_text(aes(label = paste0(Freq, "\n(", scales::percent(Freq_percentage, accuracy = 0.1), ")")),
              color = "black", size = 5) +
    scale_fill_gradientn(colors = palette_colors, limits = c(0, 1), labels = scales::percent, name = NULL) +
    scale_x_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    scale_y_discrete(labels = c("0" = "non-FHdRCC", "1" = "FHdRCC")) +
    labs(title = dataset_name, x = "Predicted", y = "Actual") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 11, color = "black"),
      panel.grid = element_blank(), 
      legend.position = "right"
    ) +
    coord_fixed()
  
  return(p)
}

# 5. Generation and Export
# ------------------------------------------------------------------------------
plot_list <- list()

for (ds in datasets) {
  sub_data <- subset(patient_data, Dataset == ds)
  cm_obj <- confusionMatrix(sub_data$pred_labels, sub_data$true_labels)
  cm_df <- as.data.frame(cm_obj$table)
  
  current_palette <- color_palettes[[ds]] %||% c("#FFFFFF", "#CCCCCC", "#000000")
  plot_list[[ds]] <- plot_cm_square(cm_df, ds, current_palette)
}

# Arrange and save
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

combined_plot <- grid.arrange(grobs = plot_list, ncol = 3)

ggsave(filename = file.path(results_dir, "Figure_S6_Patient_Confusion_Matrices.pdf"), 
       plot = combined_plot, width = 15, height = 5, device = "pdf")

message("Patient-level matrices saved to /results/.")


# ==============================================================================
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
# Script Name: Figure_S10_S12_Subgroup_ROC.R
# Description: Generate ROC curves for subgroups in Figure S10 and S12.
#              - Curve Shape: Based on 'pred_probs' (Continuous).
#              - Statistics: AUROC & 95% CI based on 'pred_labels' (Binary).
#              - Layout: Automated grid arrangement based on group count.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(gridExtra)
library(here)

# Load data using relative path
data_path <- here("data", "Figure S10 S12.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

raw_data <- read.csv(data_path)

# Ensure numeric format for columns
raw_data$true_labels <- as.numeric(as.character(raw_data$true_labels))
raw_data$pred_labels <- as.numeric(as.character(raw_data$pred_labels))
raw_data$pred_probs  <- as.numeric(as.character(raw_data$pred_probs))

groups <- unique(raw_data$GROUP)

# 2. Plot Generation Loop
# ------------------------------------------------------------------------------
plot_list <- list()
line_color <- "#1C93B5" 

for (grp in groups) {
  sub_data <- subset(raw_data, GROUP == grp)
  
  # --- Part A: Statistical Computation (Binary Labels) ---
  # AUC calculated using binary prediction labels as requested
  roc_stat_obj <- roc(sub_data$true_labels, sub_data$pred_labels, quiet = TRUE)
  
  ci_result <- tryCatch({
    ci.auc(roc_stat_obj, method = "delong")
  }, error = function(e) {
    return(c(NA, NA, NA))
  })
  
  # Format annotation text
  if (all(is.na(ci_result))) {
    label_str <- "AUC: NA"
  } else {
    label_str <- paste0("AUC: ", sprintf("%.3f", ci_result[2]), "\n",
                        "95% CI: ", sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3]))
  }
  
  # Define text position df
  anno_df <- data.frame(x = 0.5, y = 0.15, label = label_str)
  
  # --- Part B: Coordinate Extraction (Continuous Probabilities) ---
  # ROC curve shape based on predicted probabilities
  roc_curve_obj <- roc(sub_data$true_labels, sub_data$pred_probs, quiet = TRUE)
  
  ds_coords <- data.frame(
    spec = roc_curve_obj$specificities, 
    tpr  = roc_curve_obj$sensitivities,
    Group = grp
  )
  
  # --- Part C: Individual Plot Creation ---
  p <- ggplot(ds_coords, aes(x = spec, y = tpr)) +
    geom_path(color = line_color, linewidth = 0.8) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey60") +
    
    # Statistical Annotation
    annotate("text", x = 0.5, y = 0.15, label = label_str, 
             color = "black", size = 3.5, hjust = 0, lineheight = 1.2) +
    
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0.02, 0)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0.02, 0)) +
    
    labs(title = grp) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
    ) +
    coord_fixed()
  
  plot_list[[grp]] <- p
}

# 3. Grid Assembly and Export
# ------------------------------------------------------------------------------
# Calculate grid dimensions automatically
n_plots <- length(plot_list)
ncol_grid <- ceiling(sqrt(n_plots))
nrow_grid <- ceiling(n_plots / ncol_grid)

combined_plot <- grid.arrange(grobs = plot_list, ncol = ncol_grid)

# Define output path
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

ggsave(filename = file.path(results_dir, "Figure_S10_S12_ROC.pdf"), 
       plot = combined_plot, 
       width = ncol_grid * 3.5, 
       height = nrow_grid * 3.5, 
       limitsize = FALSE, device = "pdf")

message("Subgroup ROC analysis complete. Plot saved to /results/.")


# ==============================================================================
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
# Script Name: Figure_S16_NRI_Validation.R
# Description: Analysis of Net Reclassification Index (NRI) and Confusion Matrices
#              to evaluate model-assisted diagnostic improvement (Figure S16).
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(here)

# Load Data using relative path
data_path <- here("data", "Figure S16.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

demo_data <- read.table(data_path, sep = ",", header = TRUE)
true_labels <- demo_data$true_labels

# Initialize storage for reader-specific analysis
all_plots <- list()

# 2. NRI Analysis and Visualization Loop (Readers 1-7)
# ------------------------------------------------------------------------------
for (i in 1:7) {
  
  message(paste("Processing NRI analysis for Reader", i, "..."))
  
  # Extract predictions: Without AI vs. With AI
  without_ai <- demo_data[[paste0("Human", i, "_Without_WCHFHPM")]]
  with_ai    <- demo_data[[paste0("Human", i, "_With_WCHFHPM")]]
  
  # --- Part A: Confusion Matrices (Sub-group Analysis) ---
  
  # 1. Confusion Matrix for FHdRCC cases (True Positives)
  idx_1 <- which(true_labels == 1)
  cm_1  <- table(First = factor(without_ai[idx_1], levels=0:1), 
                 Second = factor(with_ai[idx_1], levels=0:1))
  
  # 2. Confusion Matrix for non-FHdRCC cases (True Negatives)
  idx_0 <- which(true_labels == 0)
  cm_0  <- table(First = factor(without_ai[idx_0], levels=0:1), 
                 Second = factor(with_ai[idx_0], levels=0:1))
  
  # --- Part B: Net Reclassification Index (NRI) Calculation ---
  
  # Counts of reclassified cases
  improved_tp <- sum(without_ai == 0 & with_ai == 1 & true_labels == 1) # FN -> TP
  worsened_fn <- sum(without_ai == 1 & with_ai == 0 & true_labels == 1) # TP -> FN
  
  improved_tn <- sum(without_ai == 1 & with_ai == 0 & true_labels == 0) # FP -> TN
  worsened_fp <- sum(without_ai == 0 & with_ai == 1 & true_labels == 0) # TN -> FP
  
  # Reclassification probabilities
  p_event_up    <- (improved_tp - worsened_fn) / sum(true_labels == 1)
  p_nonevent_up <- (improved_tn - worsened_fp) / sum(true_labels == 0)
  
  nri_total <- p_event_up + p_nonevent_up
  
  # --- Part C: Visualization ---
  
  # Prepare data frames for heatmaps
  cm1_df <- as.data.frame(as.table(cm_1))
  levels(cm1_df$First) <- levels(cm1_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  cm0_df <- as.data.frame(as.table(cm_0))
  levels(cm0_df$First) <- levels(cm0_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  # Plot 1: FHdRCC Confusion Matrix
  p1 <- ggplot(cm1_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    labs(title = "Ground truth: FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 2: non-FHdRCC Confusion Matrix
  p2 <- ggplot(cm0_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15") +
    labs(title = "Ground truth: non-FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 3: NRI Component Bar Chart
  nri_df <- data.frame(
    Category = c("FHdRCC Improvement", "non-FHdRCC Improvement", "Total NRI"),
    Value = c(p_event_up, p_nonevent_up, nri_total)
  )
  
  p3 <- ggplot(nri_df, aes(x = Category, y = Value, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    labs(title = "NRI Components", x = NULL, y = "Value") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store combined reader plot
  all_plots[[i]] <- p1 + p2 + p3 + plot_layout(ncol = 3)
}

# 3. Final Assembly and Export
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Vertically stack all reader results
final_figure <- wrap_plots(all_plots, ncol = 1)

ggsave(filename = file.path(results_dir, "Figure_S16_NRI_Validation.pdf"), 
       plot = final_figure, width = 12, height = 4 * 7, limitsize = FALSE)

message("NRI Validation analysis complete. PDF saved to /results/.")


# ==============================================================================
# Script Name: Figure_S17_NRI_Analysis.R
# Description: Analysis of Net Reclassification Index (NRI) and Confusion Matrices
#              to quantify diagnostic improvement with WCHFHPM model support.
# Author: Jinge Zhao
# ==============================================================================

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(ggplot2)
library(patchwork)
library(here)

# Use relative path for data import
data_path <- here("data", "Figure S17.csv")
if (!file.exists(data_path)) {
  stop("Data file not found! Please check /data/ directory.")
}

demo_data <- read.table(data_path, sep = ",", header = TRUE)
true_labels <- demo_data$true_labels

# Initialize list to store plot objects for all readers
all_plots <- list()

# 2. NRI Calculation and Plotting Loop (Readers 1 to 7)
# ------------------------------------------------------------------------------
for (i in 1:7) {
  
  message(paste("Analyzing performance for Reader", i, "..."))
  
  # Extract prediction vectors
  pred_without <- demo_data[[paste0("Human", i, "_Without_WCHFHPM")]]
  pred_with    <- demo_data[[paste0("Human", i, "_With_WCHFHPM")]]
  
  # --- Part A: Confusion Matrices (Comparison within Ground Truth subsets) ---
  
  # Subset 1: True FHdRCC cases
  idx_1 <- which(true_labels == 1)
  cm_1  <- table(First = factor(pred_without[idx_1], levels=0:1), 
                 Second = factor(pred_with[idx_1], levels=0:1))
  
  # Subset 2: True non-FHdRCC cases
  idx_0 <- which(true_labels == 0)
  cm_0  <- table(First = factor(pred_without[idx_0], levels=0:1), 
                 Second = factor(pred_with[idx_0], levels=0:1))
  
  # --- Part B: NRI Calculation ---
  
  # Reclassification logic
  improved_tp <- sum(pred_without == 0 & pred_with == 1 & true_labels == 1) # FN -> TP
  worsened_fn <- sum(pred_without == 1 & pred_with == 0 & true_labels == 1) # TP -> FN
  
  improved_tn <- sum(pred_without == 1 & pred_with == 0 & true_labels == 0) # FP -> TN
  worsened_fp <- sum(pred_without == 0 & pred_with == 1 & true_labels == 0) # TN -> FP
  
  # Probabilities of reclassification
  p_event_up    <- (improved_tp - worsened_fn) / sum(true_labels == 1)
  p_nonevent_up <- (improved_tn - worsened_fp) / sum(true_labels == 0)
  
  total_nri <- p_event_up + p_nonevent_up
  
  # --- Part C: Visualization ---
  
  # Prepare data for Heatmaps
  cm1_df <- as.data.frame(as.table(cm_1))
  levels(cm1_df$First) <- levels(cm1_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  cm0_df <- as.data.frame(as.table(cm_0))
  levels(cm0_df$First) <- levels(cm0_df$Second) <- c('non-FHdRCC', 'FHdRCC')
  
  # Plot 1: FHdRCC Confusion Matrix
  p1 <- ggplot(cm1_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
    labs(title = "GT: FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 2: non-FHdRCC Confusion Matrix
  p2 <- ggplot(cm0_df, aes(x = First, y = Second, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 4) +
    scale_fill_gradient(low = "#fee0d2", high = "#a50f15") +
    labs(title = "GT: non-FHdRCC", x = "Without AI", y = "With AI") +
    theme_minimal()
  
  # Plot 3: NRI Component Analysis
  nri_df <- data.frame(
    Category = c("FHdRCC Imp.", "non-FHdRCC Imp.", "Total NRI"),
    Value = c(p_event_up, p_nonevent_up, total_nri)
  )
  
  p3 <- ggplot(nri_df, aes(x = Category, y = Value, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
    labs(title = "NRI Results", x = NULL, y = "Value") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Combine reader-specific charts
  all_plots[[i]] <- p1 + p2 + p3 + plot_layout(ncol = 3)
}

# 3. Export Results
# ------------------------------------------------------------------------------
results_dir <- here("results")
if (!dir.exists(results_dir)) dir.create(results_dir)

# Vertically stack all reader panels
final_fig <- wrap_plots(all_plots, ncol = 1)

ggsave(filename = file.path(results_dir, "Figure_S17_NRI_Analysis.pdf"), 
       plot = final_fig, width = 12, height = 4 * 7, limitsize = FALSE)

message("Analysis complete. Combined PDF exported to /results/.")


# ==============================================================================
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