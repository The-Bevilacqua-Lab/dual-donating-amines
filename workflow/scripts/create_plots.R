library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tidytext)
library(ggh4x)
library(ggpattern)

#### H-BONDING MEASUREMENT PLOTS ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE, na.strings = "NaN", comment.char = "#")

# Extract the rows from the combined data frame relevant for these plots.
pairs_df <- combined_df %>% filter(geom == 1 & don_resn %in% c("A", "C", "G")) %>% 
  distinct(don_index, acc_index, eq_class_member, .keep_all = TRUE)

# Find the maximum fill value corresponding to the general region below the
# h_dist_max cutoff specified in the config file.
h_bond_region <- pairs_df[(pairs_df$h_acc_distance <= snakemake@config[["h_dist_max"]]),] %>%
  ggplot(aes(x = h_angle, y = h_acc_distance)) + geom_bin_2d(binwidth = c(120/100, 2.0/100))
max_value <- max(ggplot_build(h_bond_region)$data[[1]]$value)

# Create the plot.
pairs <- pairs_df %>% ggplot(aes(x = h_angle, y = h_acc_distance)) +
  geom_bin_2d(binwidth = c(120/100, 2.0/100)) + 
  geom_segment(x=140, y=1.0, xend=140, yend=2.5, linewidth=0.4, linetype=2, colour="Red") +
  geom_segment(x=140, y=2.5, xend=180, yend=2.5, linewidth=0.4, linetype=2, colour="Red") +
  geom_segment(x=180, y=1.0, xend=180, yend=2.5, linewidth=0.4, linetype=2, colour="Red") +
  geom_segment(x=140, y=1.0, xend=180, yend=1.0, linewidth=0.4, linetype=2, colour="Red") +
  scale_fill_viridis(limits = c(1, max_value), name = "Count") +
  xlab("Angle (\ub0)") +
  ylab("Distance (\uc5)") +
  coord_fixed(ratio = 120/2.0, xlim = c(60, 180), ylim = c(1.0, 3.0)) +
  scale_x_continuous(breaks = seq(60, 180, 30)) +
  scale_y_continuous(breaks = seq(1.0, 3.0, 0.4)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.key.width = unit(0.25, "in"), legend.key.height = unit(0.125, "in"),
        legend.title = element_text(vjust = 0.9),
        axis.title.x = element_text(hjust = 0.29, color = rgb(0/255, 158/255, 115/255),
                                    margin = margin(t = 0.25, r = 0, b = 0.125, l = 0, unit = "in")),
        axis.title.y = element_text(hjust = 0.26, color = rgb(213/255, 94/255, 0/255),
                                    margin = margin(t = 0, r = 0.25, b = 0, l = 0.125, unit = "in")))

# Write the plot.
ggsave(snakemake@output[["pairs"]], plot = pairs, width = 3.25, height = 3.75, units = "in", scale = 1)

#### PSEUDOTORSION PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
pt_df <- combined_df %>%
  filter(don_resn %in% c("A", "C", "G") & eta != -360 & theta != -360) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
pt_df[pt_df$type == 0, "type"] <- "No"
pt_df[pt_df$type == 1, "type"] <- "Single"
pt_df[pt_df$type == 2, "type"] <- "Dual"
pt_df$type <- factor(pt_df$type, levels = c("No", "Single", "Dual"))

# Translate the pseudotorsions to range from 0 to 360 degrees.
pt_df[which(pt_df$eta >= 0), "eta_translated"] <- pt_df[which(pt_df$eta >= 0), "eta"]
pt_df[which(pt_df$eta < 0), "eta_translated"] <- pt_df[which(pt_df$eta < 0), "eta"] + 360
pt_df[which(pt_df$theta >= 0), "theta_translated"] <- pt_df[which(pt_df$theta >= 0), "theta"]
pt_df[which(pt_df$theta < 0), "theta_translated"] <- pt_df[which(pt_df$theta < 0), "theta"] + 360

# Set the width of the square bin.
bin_w <- 360/100

# Set the boundaries of the central region.
eta_min <- 90+bin_w*16
eta_max <- 180+bin_w*2
theta_min <- 90+bin_w*21
theta_max <- 180+bin_w*24

# Prepare a data frame with pseudotorsions outside the central region.
pt_df_out <- pt_df %>%
  filter(((eta_translated < eta_min | eta_translated > eta_max) &
            (theta_translated > theta_min | theta_translated < theta_max)) |
           ((theta_translated < theta_min | theta_translated > theta_max) &
              (eta_translated > eta_min | eta_translated < eta_max)))

# Create the plot with all eta theta values and a red-dashed box outlining the central region.
pt_no_all <- pt_df[which(pt_df$type == "No"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.1, linetype=2, colour="Red", fill=NA) +
  scale_fill_viridis(name = "Count") +
  ggtitle("No") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_single_all <- pt_df[which(pt_df$type == "Single"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.1, linetype=2, colour="Red", fill=NA) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Single") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_dual_all <- pt_df[which(pt_df$type == "Dual"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.1, linetype=2, colour="Red", fill=NA) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Dual") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_plot_all <- grid.arrange(pt_no_all, pt_single_all, pt_dual_all, nrow = 3)

# Create the plot with eta theta values outside the central region.
pt_no_out <- pt_df_out[which(pt_df_out$type == "No"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("No") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_single_out <- pt_df_out[which(pt_df_out$type == "Single"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Single") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_dual_out <- pt_df_out[which(pt_df_out$type == "Dual"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Dual") +
  xlab("Eta (\ub0)") +
  ylab("Theta (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_plot_out <- grid.arrange(pt_no_out, pt_single_out, pt_dual_out, nrow = 3)

# Write the plot.
ggsave(snakemake@output[["pseudotorsion_all"]], plot = pt_plot_all, width = 3.25, height = 6, units = "in", scale = 1)
ggsave(snakemake@output[["pseudotorsion_out"]], plot = pt_plot_out, width = 3.25, height = 6, units = "in", scale = 1)

#### HEAVY ATOM DENSITY PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
density_all_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
density_all_df[density_all_df$type == 0, "type"] <- "No"
density_all_df[density_all_df$type == 1, "type"] <- "Single"
density_all_df[density_all_df$type == 2, "type"] <- "Dual"
density_all_df$type <- factor(density_all_df$type, levels = c("No", "Single", "Dual"))

# Create a column that contains the donor atom name along with the donor residue name.
density_all_df <- density_all_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Add density columns.
volume_1 <- (4/3)*pi*snakemake@config[["count_dist_1"]]**3
volume_2 <- (4/3)*pi*snakemake@config[["count_dist_2"]]**3
volume_diff <- volume_2 - volume_1
density_df$density_1 <- density_df$count_1/volume_1
density_df$density_2 <- (density_df$count_2-density_df$count_1)/volume_diff

# Pivot the density_1 and density_2 columns to place the values in a single column.
density_all_pivot_df <- density_all_df %>%
  pivot_longer(cols = c(density_1, density_2), names_to = "density_type", values_to = "density_values")

# Change the density_type naming.
density_all_pivot_df[density_all_pivot_df$density_type == "density_1", "density_type"] <- "ROI 1"
density_all_pivot_df[density_all_pivot_df$density_type == "density_2", "density_type"] <- "ROI 2"

# Calculate sample sizes from the non-filtered data set.
density_summary_all <- density_all_pivot_df %>%
  summarise(n = n(), den_med = median(density_values), .by = c(type, don_label, density_type))

# Color Palette: black, orange, sky blue, bluish green, yellow, blue, vermilion, reddish purple
# Colors from DOI 10.1038/nmeth.1618
color_palette <- c(rgb(0/255, 0/255, 0/255), rgb(230/255, 159/255, 0/255), rgb(86/255, 180/255, 233/255),
                   rgb(0/255, 158/255, 115/255), rgb(240/255, 228/255, 66/255), rgb(0/255, 114/255, 178/255),
                   rgb(213/255, 94/255, 0/255), rgb(204/255, 121/255, 167/255))

# Specify the variables for the plots.
ylim <- c(-80, 200)
ratio <- (3/diff(ylim))*2

# Create the box plots with all data (outliers included).
density_all_plot <- density_all_pivot_df %>% ggplot(aes(x=type, y=density_values*1000, fill=density_type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, show.legend = FALSE) +
  geom_text(data = density_summary_all,
            aes(x=type, y=-45, label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 8, size.unit = "pt", hjust = 0.5, vjust = 0.5, angle = 60) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  scale_fill_manual(values = c(color_palette[3], color_palette[6])) +
  xlab("Type of Amine H-Bonding") +
  ylab("Density (heavy-atoms/nm\ub3)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_nested_wrap(vars(density_type, don_label), nrow = 1, scales = "fixed")

# Write the plots.
ggsave(snakemake@output[["density_all"]], plot = density_all_plot, width = 6, height = 3, units = "in", scale = 1)

# Calculate quantile values for both densities.
density_quantiles <- density_all_df %>% summarise(iqr_1 = IQR(density_1),
                                                  q1_1 = quantile(density_1, 0.25),
                                                  q3_1 = quantile(density_1, 0.75),
                                                  iqr_2 = IQR(density_2),
                                                  q1_2 = quantile(density_2, 0.25),
                                                  q3_2 = quantile(density_2, 0.75),
                                                  .by = c(type, don_label))

# Calculate some stats on the IQR for density 1.
avg_iqr_1_n <- sum(density_quantiles[density_quantiles$type == "No",]$iqr_1)/3 # 0.0154675
avg_iqr_1_s <- sum(density_quantiles[density_quantiles$type == "Single",]$iqr_1)/3 # 0.009242776
avg_iqr_1_d <- sum(density_quantiles[density_quantiles$type == "Dual",]$iqr_1)/3 # 0.007922379
n_d_1_drop <- ((avg_iqr_1_n-avg_iqr_1_d)/avg_iqr_1_n)*100 # 48.78049 %
s_d_1_drop <- ((avg_iqr_1_s-avg_iqr_1_d)/avg_iqr_1_s)*100 # 14.28571 %

# Calculate some stats on the IQR for density 2.
avg_iqr_2_n <- sum(density_quantiles[density_quantiles$type == "No",]$iqr_2)/3 # 0.01392768
avg_iqr_2_s <- sum(density_quantiles[density_quantiles$type == "Single",]$iqr_2)/3 # 0.01164105
avg_iqr_2_d <- sum(density_quantiles[density_quantiles$type == "Dual",]$iqr_2)/3 # 0.008782754
n_d_2_drop <- ((avg_iqr_2_n-avg_iqr_2_d)/avg_iqr_2_n)*100 # 36.9403 %
s_d_2_drop <- ((avg_iqr_2_s-avg_iqr_2_d)/avg_iqr_2_s)*100 # 24.55357 %

# Create a data frame with outliers for both densities filtered out.
density_filtered_df <- merge(density_all_df, density_quantiles) %>%
  filter(density_1 >= q1_1 - 2 * iqr_1 & density_1 <= q3_1 + 2 * iqr_1) %>%
  filter(density_2 >= q1_2 - 2 * iqr_2 & density_2 <= q3_2 + 2 * iqr_2)

# Pivot the density_1 and density_2 columns to place the values in a single column.
density_filtered_pivot_df <- density_filtered_df %>%
  pivot_longer(cols = c(density_1, density_2), names_to = "density_type", values_to = "density_values")

# Change the density_type naming.
density_filtered_pivot_df[density_filtered_pivot_df$density_type == "density_1", "density_type"] <- "ROI 1"
density_filtered_pivot_df[density_filtered_pivot_df$density_type == "density_2", "density_type"] <- "ROI 2"

# Calculate sample sizes from the filtered data set.
density_summary_filtered <- density_filtered_pivot_df %>%
  summarise(n = n(), den_med = median(density_values), .by = c(type, don_label, density_type))

# Create the labels for sample sizes. Labels for ROI 2 will be blank to save room on the plots. The sample sizes for
# ROI 2 are the same as those for ROI 1.
density_summary_filtered <- density_summary_filtered %>%
  mutate(sample_label = paste("n = ", prettyNum(n, big.mark = ","), sep = ""))
density_summary_filtered[density_summary_filtered$density_type == "ROI 2", "sample_label"] <- ""

# Specify the variables for the plots.
ylim <- c(-15, 80)
ratio <- (3/diff(ylim))*2

# Create the box plots with outliers filtered out.
density_filtered_plot <- density_filtered_pivot_df %>% ggplot(aes(x=type, y=density_values*1000, fill=density_type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, show.legend = FALSE) +
  geom_text(data = density_summary_filtered,
            aes(x=type, y=-4, label = sample_label),
            size = 8, size.unit = "pt", hjust = 0.5, vjust = 0.5, angle = 60) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  scale_fill_manual(values = c(color_palette[3], color_palette[6])) +
  xlab("Type of Amine H-Bonding") +
  ylab("Density (heavy-atoms/nm\ub3)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_nested_wrap(vars(density_type, don_label), nrow = 1, scales = "fixed")

# Write the plots.
ggsave(snakemake@output[["density_filtered"]], plot = density_filtered_plot, width = 6, height = 3, units = "in", scale = 1)

# Perform Mann-Whitney tests for density 2 with outliers removed.

# I do not think the distributions being compared are fully independent from one another, an assumption of the
# Mann-Whitney test. I may want to consider this further before confidently relying on these p-values.

# Adenine
adenine <- density_filtered_df[(density_filtered_df$don_resn == "A"),]
a_no <- adenine[adenine$type == "No", "density_2"]
a_single <- adenine[adenine$type == "Single", "density_2"]
a_dual <- adenine[adenine$type == "Dual", "density_2"]
# A - No to A - Single
cat("\nA - No to A - Single Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(a_no, a_single)) # W = 352982486, p-value < 2.2e-16
# A - Single to A - Dual
cat("\nA - Single to A - Dual Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(a_single, a_dual)) # W = 47115675, p-value < 2.2e-16

# Cytosine
cytosine <- density_filtered_df[(density_filtered_df$don_resn == "C"),]
c_no <- cytosine[cytosine$type == "No", "density_2"]
c_single <- cytosine[cytosine$type == "Single", "density_2"]
c_dual <- cytosine[cytosine$type == "Dual", "density_2"]
# C - No to C - Single
cat("\nC - No to C - Single Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(c_no, c_single)) # W = 218199061, p-value < 2.2e-16
# C - Single to C - Dual
cat("\nC - Single to C - Dual Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(c_single, c_dual)) # W = 35739526, p-value < 2.2e-16

# Guanine
guanine <- density_filtered_df[(density_filtered_df$don_resn == "G"),]
g_no <- guanine[guanine$type == "No", "density_2"]
g_single <- guanine[guanine$type == "Single", "density_2"]
g_dual <- guanine[guanine$type == "Dual", "density_2"]
# G - No to G - Single
cat("\nG - No to G - Single Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(g_no, g_single)) # W = 400955757, p-value = 2.059e-06
# G - Single to G - Dual
cat("\nG - Single to G - Dual Mann-Whitney test for density 2 with outliers removed\n\n")
print(wilcox.test(g_single, g_dual)) # W = 124882421, p-value < 2.2e-16

#### SASA PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
sasa_all_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
sasa_all_df[sasa_all_df$type == 0, "type"] <- "No"
sasa_all_df[sasa_all_df$type == 1, "type"] <- "Single"
sasa_all_df[sasa_all_df$type == 2, "type"] <- "Dual"
sasa_all_df$type <- factor(sasa_all_df$type, levels = c("No", "Single", "Dual"))

# Create a column that contains the donor atom name along with the donor residue name.
sasa_all_df <- sasa_all_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Calculate sample sizes from the non-filtered data set.
sasa_summary_all <- sasa_all_df %>%
  summarise(n = n(), sasa_med = median(SASA), .by = c(type, don_label))

# Check the median values for C single donating amines and G no donating amines.
c_single <- sasa_all_df %>% filter(type == "Single" & don_resn == "C")
print(median(c_single$SASA)) # 8.748707
g_no <- sasa_all_df %>% filter(type == "No" & don_resn == "G")
print(median(g_no$SASA)) # 8.748707

# Specify the variables for the plots.
ylim <- c(-20, 60)
ratio <- (3/diff(ylim))*2

# Create the box plots with all data (outliers included).
sasa_all_box_plot <- sasa_all_df %>% ggplot(aes(x=type, y=SASA)) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = sasa_summary_all, aes(x=type, y=-12, label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 8, size.unit = "pt", hjust = 0.5, vjust = 0.5, angle = 60) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  xlab("Type of Amine H-Bonding") +
  ylab("SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_label, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["sasa_all_box"]], plot = sasa_all_box_plot, width = 3.5, height = 4, units = "in", scale = 1)

# Specify the variables for the plots.
ylim <- c(0, 15)
ratio <- (3/diff(ylim))*2

# Create the column plot with all data (outliers included).
sasa_all_col_plot <- sasa_summary_all %>% ggplot(aes(x=type, y=sasa_med)) +
  geom_col(width = 0.8, linewidth = 0.3, color = "black", fill = "grey", show.legend = FALSE) +
  geom_text(aes(x=type, y=sasa_med+max(sasa_med)*0.05, label=round(sasa_med, digits = 1)),
            size = 10, size.unit = "pt", vjust = 0, inherit.aes = FALSE) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  xlab("Type of Amine H-Bonding") +
  ylab("Median SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap( ~ don_label, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["sasa_all_col"]], plot = sasa_all_col_plot, width = 3.5, height = 4, units = "in", scale = 1)

# Calculate quantile values for both densities.
sasa_quantiles <- sasa_all_df %>% summarise(iqr = IQR(SASA),
                                            q1 = quantile(SASA, 0.25),
                                            q3 = quantile(SASA, 0.75),
                                            .by = c(type, don_label))

# Calculate some stats on the IQRs.
avg_iqr_n <- sum(sasa_quantiles[sasa_quantiles$type == "No",]$iqr)/3 # 18.591
avg_iqr_s <- sum(sasa_quantiles[sasa_quantiles$type == "Single",]$iqr)/3 # 9.842296
avg_iqr_d <- sum(sasa_quantiles[sasa_quantiles$type == "Dual",]$iqr)/3 # 0.3645295
n_d_drop <- ((avg_iqr_n-avg_iqr_d)/avg_iqr_n)*100 # 98.03922 %
s_d_drop <- ((avg_iqr_s-avg_iqr_d)/avg_iqr_s)*100 # 96.2963 %

# Create a data frame with outliers for both densities filtered out.
sasa_filtered_df <- merge(sasa_all_df, sasa_quantiles) %>%
  filter(SASA >= q1 - 2 * iqr & SASA <= q3 + 2 * iqr)

# Calculate sample sizes from the filtered data set.
sasa_summary_filtered <- sasa_filtered_df %>%
  summarise(n = n(), sasa_med = median(SASA), .by = c(type, don_label))

# Specify the variables for the plots.
ylim <- c(-8.5, 60)
ratio <- (3/diff(ylim))*1.5

# Create the box plots with outliers filtered out.
sasa_filtered_box_plot <- sasa_filtered_df %>% ggplot(aes(x=type, y=SASA)) +
  geom_violin(fill = "grey", show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = sasa_summary_filtered, aes(x=type, y=-4,
                                              label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 10, size.unit = "pt", vjust = 1, angle = 30) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  xlab("Type of Amine H-Bonding") +
  ylab("SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_label, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["sasa_filtered_box"]], plot = sasa_filtered_box_plot, width = 6.5, height = 9, units = "in",
       scale = 1)

#### DONOR IDENTITY PLOT ####

# Extract the rows from the combined data frame relevant for these plots.
don_id_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
don_id_df[don_id_df$type == 0, "type"] <- "No"
don_id_df[don_id_df$type == 1, "type"] <- "Single"
don_id_df[don_id_df$type == 2, "type"] <- "Dual"
don_id_df$type <- factor(don_id_df$type, levels = c("No", "Single", "Dual"))

# Create a column that contains the donor atom name along with the donor residue name.
don_id_df <- don_id_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Calculate samples sizes and merge into dataframe.
don_id_df <- merge(don_id_df, summarise(don_id_df, n_resn = n(), .by = c(don_label)))
don_id_df <- merge(don_id_df, summarise(don_id_df, n_don_type = n(), .by = c(don_label, type)))

# Only keep one row for each donor residue and donor type combination.
don_id_df <- don_id_df %>% distinct(don_label, type, .keep_all = TRUE)

# Add a column that specifies the percent occurrence for each category.
don_id_df <-don_id_df %>% mutate(occurance = n_don_type/sum(n_don_type)*100)

# Create the plot.
don_id_plot <- don_id_df %>% ggplot(aes(x=type, y=occurance)) +
  geom_col(width = 0.8, linewidth = 0.3, color = "black", fill = "grey", show.legend = FALSE) +
  geom_text(aes(x=type, y=occurance+max(occurance)*0.05, label=paste(round(occurance, digits = 1), "%", sep = "")),
            vjust = 0, inherit.aes = FALSE) +
  xlab("Type of Amine H-Bonding") +
  ylab("Occurrence (%)") +
  scale_y_continuous(limits = c(0, 27)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1) +
  facet_wrap( ~ don_label, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["don_id"]], plot = don_id_plot, width = 6, height = 5, units = "in", scale = 1)

#### ACCEPTOR PAIR IDENTITY PLOT ####

# Extract the rows from the combined data frame relevant for these plots.
acc_pair_id_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G") & type == 2) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Create a column that contains the donor atom name along with the donor residue name.
acc_pair_id_df <- acc_pair_id_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Create a data frame with just acceptor pairs that belong to the same residue.
acc_pair_id_same_df <- acc_pair_id_df %>% filter(same_resi == "True")

# Calculate samples sizes and merge into the data frames.
acc_pair_id_df <- merge(acc_pair_id_df, summarise(acc_pair_id_df, n_resn_pair = n(),
                                                  .by = c(don_label, acc_pair_combined)))
acc_pair_id_same_df <- merge(acc_pair_id_same_df, summarise(acc_pair_id_same_df, n_resn_pair_same = n(),
                                                            .by = c(don_label, acc_pair_combined)))

# Only keep one row for each donor resn and acceptor pair combination.
acc_pair_id_df <- acc_pair_id_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)
acc_pair_id_same_df <- acc_pair_id_same_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)

# Add a column specifying the sample sizes for acceptor pairs that belong to the same residue.
acc_pair_id_df <- merge(acc_pair_id_df, acc_pair_id_same_df[,c("don_label", "acc_pair_combined", "n_resn_pair_same")],
                        all.x = TRUE)

# Create a column that indicates the fraction of acceptor pairs that belong to the same residue. Fractions that would
# be 0 will be given a value of NA.
acc_pair_id_df <- acc_pair_id_df %>% mutate(fraction_label = NA)
acc_pair_id_df[which(!is.na(acc_pair_id_df$n_resn_pair_same)), "fraction_label"] <-
  paste(acc_pair_id_df[which(!is.na(acc_pair_id_df$n_resn_pair_same)), "n_resn_pair_same"], "/",
        acc_pair_id_df[which(!is.na(acc_pair_id_df$n_resn_pair_same)), "n_resn_pair"], sep = " ")

# Create the plot.
acc_pair_id_plot <- acc_pair_id_df %>%
  ggplot(aes(x=reorder_within(acc_pair_combined, n_resn_pair, don_label), y=n_resn_pair)) +
  geom_col(width = 0.8, color = "black", fill = "grey", show.legend = FALSE) +
  geom_col_pattern(aes(y=n_resn_pair_same),
                   width = 0.8, color = "black", fill = "grey", pattern = "stripe", pattern_fill = "black",
                   pattern_size = 0, pattern_angle = 45, pattern_spacing = 0.02) +
  geom_text(aes(y=n_resn_pair+max(n_resn_pair)*0.05, label=fraction_label),
            size = 10, size.unit = "pt", hjust = 0, vjust = 0, angle = 45) +
  xlab("Pair of Acceptor Atoms") +
  ylab("Count") +
  scale_x_reordered(limits = function(x) rev(x)[1:10]) +
  scale_y_continuous(limits = c(0, 2100)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap( ~ don_label, nrow = 1, scales = "free_x")

# Write the plots.
ggsave(snakemake@output[["acc_pair_id"]], plot = acc_pair_id_plot, width = 6, height = 5, units = "in", scale = 1)

### CHI PLOT ###

# Extract the rows from the combined data frame relevant for these plots.
chi_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
chi_df[chi_df$type == 0, "type"] <- "No"
chi_df[chi_df$type == 1, "type"] <- "Single"
chi_df[chi_df$type == 2, "type"] <- "Dual"
chi_df$type <- factor(chi_df$type, levels = c("No", "Single", "Dual"))

# Adjust the chi dihedrals such that values of -180 are instead 180.
chi_df <- chi_df %>% mutate(chi_adjusted = chi)
chi_df[which(chi_df$chi_adjusted == -180), "chi_adjusted"] <- 180

# Write out the full names of the nucleobases.
chi_df[which(chi_df$don_resn == "A"), "don_resn"] <- "Adenine"
chi_df[which(chi_df$don_resn == "C"), "don_resn"] <- "Cytosine"
chi_df[which(chi_df$don_resn == "G"), "don_resn"] <- "Guanine"

# Split the chi_df into groups based on don_resn and type.
chi_split <- split(chi_df, ~ don_resn + type)

# Create a data frame with chi data binned and grouped by don_resn and type.
chi_bins <- data.frame(don_resn = character(), type = character(), counts = integer(), density = double(), 
                       mids = integer())
for (split_idx in seq_along(chi_split)) {
  group_bin <- hist(chi_split[[split_idx]]$chi_adjusted, breaks = seq(-180, 180, 10), plot = FALSE)
  for (hist_idx in seq_along(group_bin$mids)) {
    chi_bins <- rbind(chi_bins, data.frame(
      don_resn = unlist(strsplit(names(chi_split)[split_idx], "[.]"))[1],
      type = unlist(strsplit(names(chi_split)[split_idx], "[.]"))[2],
      counts = group_bin$counts[hist_idx],
      density = group_bin$density[hist_idx],
      mids = group_bin$mids[hist_idx]))
  }
}

# Factor the type column in chi_bins.
chi_bins$type <- factor(chi_bins$type, levels = c("No", "Single", "Dual"))

# Create the plots of a segment of the 360 range.
chi_plot_partial <- chi_bins %>% ggplot(aes(x=mids, y=density)) +
  geom_col(fill = "grey", show.legend = FALSE) +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi/8, end = pi*3/4) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) +
  scale_y_continuous(limits = c(0, 0.0018)) +
  xlab("\u03c7 (\ub0)") +
  ylab(element_blank()) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.3, "in"),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_nested_wrap(vars(don_resn, type), nrow = 3, scales = "fixed")

# Create the plots of a segment of the 360 range with the y-axis displayed.
chi_plot_partial_y <- chi_bins %>% ggplot(aes(x=mids, y=density)) +
  geom_col(fill = "grey", show.legend = FALSE) +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi/8, end = pi*3/4) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) +
  scale_y_continuous(limits = c(0, 0.0018)) +
  xlab("\u03c7 (\ub0)") +
  ylab(element_blank()) +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.3, "in"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_nested_wrap(vars(don_resn, type), nrow = 3, scales = "fixed")

# Create a plot showing all data combined.
chi_plot_combined <- chi_df %>% ggplot(aes(x=chi_adjusted)) +
  geom_histogram(aes(y=after_stat(count)), binwidth = 10, center = 5, show.legend = FALSE) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 1, xend = -20, yend = 10^5) +
  annotate("segment", color = "red", linetype = "dashed", x = 130, y = 1, xend = 130, yend = 10^5) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 1, xend = 130, yend = 1) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 10^5, xend = 130, yend = 10^5) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 45)) +
  scale_y_continuous(transform = "log10", name = "Count") +
  xlab("\u03c7 (\ub0)") +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 0.25)

# Write the plots.
ggsave(snakemake@output[["chi_partial"]], plot = chi_plot_partial, width = 3, height = 9, units = "in", scale = 1)
ggsave(snakemake@output[["chi_partial_y"]], plot = chi_plot_partial_y, width = 3, height = 9, units = "in", scale = 1)
ggsave(snakemake@output[["chi_combined"]], plot = chi_plot_combined, width = 6.5, height = 9, units = "in", scale = 1)

# Remove the default Rplots.pdf generated.
file.remove('Rplots.pdf')
