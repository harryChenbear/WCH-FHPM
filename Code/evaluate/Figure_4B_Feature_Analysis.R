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
