library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tidytext)
library(ggh4x)
library(ggpattern)
library(shadowtext)
library(scales)
library(cowplot)
library(circular)

#### SET INITIAL VARIABLES ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE, na.strings = "NaN", comment.char = "#")

# Color Palette: black, orange, sky blue, bluish green, yellow, blue, vermilion, reddish purple
# Colors from DOI 10.1038/nmeth.1618
color_palette <- c(rgb(0/255, 0/255, 0/255), rgb(230/255, 159/255, 0/255), rgb(86/255, 180/255, 233/255),
                   rgb(0/255, 158/255, 115/255), rgb(240/255, 228/255, 66/255), rgb(0/255, 114/255, 178/255),
                   rgb(213/255, 94/255, 0/255), rgb(204/255, 121/255, 167/255))

#### H-BONDING MEASUREMENT PLOTS ####

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
  geom_segment(x=140, y=1.0, xend=140, yend=2.5, linewidth=0.4, linetype=2, colour="red") +
  geom_segment(x=140, y=2.5, xend=180, yend=2.5, linewidth=0.4, linetype=2, colour="red") +
  geom_segment(x=180, y=1.0, xend=180, yend=2.5, linewidth=0.4, linetype=2, colour="red") +
  geom_segment(x=140, y=1.0, xend=180, yend=1.0, linewidth=0.4, linetype=2, colour="red") +
  scale_fill_viridis(limits = c(1, max_value), name = "Count") +
  xlab("Angle (\ub0)") +
  ylab("Distance (\uc5)") +
  coord_fixed(ratio = 120/2.0, xlim = c(60, 180), ylim = c(1.0, 3.0)) +
  scale_x_continuous(breaks = seq(60, 180, 30)) +
  scale_y_continuous(breaks = seq(1.0, 3.0, 0.4)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.key.width = unit(0.25, "in"), legend.key.height = unit(0.125, "in"),
        legend.margin = margin(t = 0.115, r = 0, b = 0, l = 0, "in"), legend.title = element_text(vjust = 0.9),
        axis.title.x = element_text(hjust = 0.29, color = rgb(0/255, 158/255, 115/255),
                                    margin = margin(t = 0.25, r = 0, b = 0.125, l = 0, unit = "in")),
        axis.title.y = element_text(hjust = 0.26, color = rgb(213/255, 94/255, 0/255),
                                    margin = margin(t = 0, r = 0.25, b = 0, l = 0.125, unit = "in")))

# Write the plot.
ggsave(snakemake@output[["pairs"]], plot = pairs, width = 3.25, height = 3.75, units = "in", scale = 1)

#### DONOR IDENTITY PLOT ####

# Extract the rows from the combined data frame relevant for these plots.
don_id_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
don_id_df[don_id_df$type == 0, "type"] <- "Non"
don_id_df[don_id_df$type == 1, "type"] <- "Single"
don_id_df[don_id_df$type == 2, "type"] <- "Dual"
don_id_df$type <- factor(don_id_df$type, levels = c("Non", "Single", "Dual"))

# Calculate samples sizes and merge into dataframe.
don_id_df <- merge(don_id_df, summarise(don_id_df, n_resn = n(), .by = c(don_resn)))
don_id_df <- merge(don_id_df, summarise(don_id_df, n_don_type = n(), .by = c(don_resn, type)))

# Only keep one row for each donor residue and donor type combination.
don_id_df <- don_id_df %>% distinct(don_resn, type, .keep_all = TRUE)

# Add a column that specifies the percent occurrence for each category.
don_id_df <-don_id_df %>% mutate(occurance = n_don_type/sum(n_don_type)*100)

# Create a column that contains the donor atom and residue names along with the number of occurrences.
don_id_df <- don_id_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", ", n = ",
                                                    prettyNum(n_resn, big.mark = ","), sep = ""))

# Create the plot.
don_id_plot <- don_id_df %>% ggplot(aes(x=type, y=occurance)) +
  geom_col(width = 0.8, linewidth = 0.3, color = "black", fill = "grey", show.legend = FALSE) +
  geom_text(aes(x=type, y=occurance+max(occurance)*0.05, label=paste(round(occurance, digits = 1), "%", sep = "")),
            size = 10, size.unit = "pt", vjust = 0, inherit.aes = FALSE) +
  geom_text(aes(x=type, y=-5, label = paste("n = ", prettyNum(n_don_type, big.mark = ","), sep = "")),
            size = 8, size.unit = "pt", hjust = 0.5, vjust = 0.5, angle = 30) +
  xlab("Type of Donating Amine") +
  ylab("Occurrence (%)") +
  scale_y_continuous(limits = c(-8, 27)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1, panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 9)) +
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

# Create a column that contains the combined acceptor pair divided by a forward slash instead of a comma.
acc_pair_id_df <- acc_pair_id_df %>%
  mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))

# Create the plot.
acc_pair_id_plot <- acc_pair_id_df %>%
  ggplot(aes(x=reorder_within(acc_pair_combined_reformat, n_resn_pair, don_label), y=n_resn_pair)) +
  geom_col(width = 0.8, color = "black", fill = "grey", show.legend = FALSE) +
  geom_col_pattern(aes(y=n_resn_pair_same),
                   width = 0.8, color = "black", fill = "grey", pattern = "stripe", pattern_fill = "black",
                   pattern_size = 0, pattern_angle = 45, pattern_spacing = 0.02) +
  geom_text(aes(y=n_resn_pair+max(n_resn_pair)*0.05, label=fraction_label),
            size = 8, size.unit = "pt", hjust = 0, vjust = 0, angle = 45) +
  xlab("Pair of Acceptor Atoms") +
  ylab("Count") +
  scale_x_reordered(limits = function(x) rev(x)[1:10]) +
  scale_y_continuous(limits = c(0, 2500)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        strip.text = element_text(size = 9)) +
  facet_wrap( ~ don_label, nrow = 1, scales = "free_x")

# Write the plots.
ggsave(snakemake@output[["acc_pair_id"]], plot = acc_pair_id_plot, width = 6, height = 5, units = "in", scale = 1)

# Write key data to a csv file.
acc_pair_id_df_csv <- acc_pair_id_df[c("don_label", "acc_pair_combined_reformat", "n_resn_pair", "n_resn_pair_same")]
write.csv(acc_pair_id_df_csv, file = snakemake@output[["acc_pair_id_csv"]], row.names = FALSE)

#### SASA PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
sasa_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
sasa_df[sasa_df$type == 0, "type"] <- "Non"
sasa_df[sasa_df$type == 1, "type"] <- "Single"
sasa_df[sasa_df$type == 2, "type"] <- "Dual"
sasa_df$type <- factor(sasa_df$type, levels = c("Non", "Single", "Dual"))

# Create a column that contains the donor atom name along with the donor residue name.
sasa_df <- sasa_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Create box plots.
sasa_box_plot <- sasa_df %>% ggplot(aes(x=type, y=SASA)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.stroke = 0, outlier.size = 1) +
  xlab("Type of Donating Amine") +
  ylab("SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank()) +
  facet_wrap( ~ don_label, nrow = 1)

# Create jitter plots of just the dual-donating amines.
sasa_dual_df <- sasa_df[sasa_df["type"] == "Dual",]
possible_values <- round(seq(0, ceiling(max(sasa_dual_df$SASA)), sort(unique(sasa_dual_df$SASA))[2]), digits = 1)
sasa_jitter_plot <- sasa_dual_df %>% ggplot(aes(x=don_label, y=SASA)) +
  geom_point(size = 0.5, alpha = 0.5, stroke = 0, position = position_jitter(seed = 1)) +
  scale_y_continuous(breaks = possible_values[seq(1, length(possible_values), 2)], minor_breaks = possible_values) +
  ggtitle("Dual-Donating Amines") +
  xlab("Nucleobase") +
  ylab("SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major.x = element_blank())

# Write the plots.
ggsave(snakemake@output[["sasa_box"]], plot = sasa_box_plot, width = 3.5, height = 3, units = "in", scale = 1)
ggsave(snakemake@output[["sasa_jitter"]], plot = sasa_jitter_plot, width = 2.5, height = 3, units = "in", scale = 1)

# Calculate some statistics.
sasa_stats <- sasa_df %>% summarise(n = n(),
                                    med = median(SASA),
                                    iqr = IQR(SASA),
                                    q1 = quantile(SASA, 0.25),
                                    q3 = quantile(SASA, 0.75),
                                    .by = c(type, don_label))

# Specify the variables for the column plots.
ylim <- c(0, max(sasa_stats$med)+max(sasa_stats$med)*0.125)
ratio <- (3/diff(ylim))*2

# Create column plots of SASA median values.
sasa_col_plot <- sasa_stats %>% ggplot(aes(x=type, y=med)) +
  geom_col(width = 0.8, linewidth = 0.3, color = "black", fill = "grey") +
  geom_text(aes(x=type, y=med+max(med)*0.075, label=round(med, digits = 1)),
            size = 10, size.unit = "pt", inherit.aes = FALSE) +
  coord_fixed(ratio = ratio, ylim = ylim) +
  xlab("Type of Donating Amine") +
  ylab("Median SASA (\uc5\ub2)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid.major.x = element_blank()) +
  facet_wrap( ~ don_label, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["sasa_col"]], plot = sasa_col_plot, width = 3.5, height = 3, units = "in", scale = 1)

# Write SASA statistics to a csv file.
write.csv(sasa_stats, file = snakemake@output[["sasa_stats"]], row.names = FALSE)

#### SASA WILCOXON RANK SUM TESTS ####

# Adenine p-values
a_sasa_df <- sasa_df %>% filter(don_resn == "A")
a_sasa_pvals <- data.frame(don_resn = "A",
                            n.s = wilcox.test(pull(filter(a_sasa_df, type == "Non"), SASA),
                                              pull(filter(a_sasa_df, type == "Single"), SASA))$p.value,
                            n.d = wilcox.test(pull(filter(a_sasa_df, type == "Non"), SASA),
                                              pull(filter(a_sasa_df, type == "Dual"), SASA))$p.value,
                            s.d = wilcox.test(pull(filter(a_sasa_df, type == "Single"), SASA),
                                              pull(filter(a_sasa_df, type == "Dual"), SASA))$p.value)
a_sasa_adjust <- p.adjust(as.numeric(a_sasa_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
a_sasa_pvals <-
  cbind(a_sasa_pvals, n.s_adjust = a_sasa_adjust[1], n.d_adjust = a_sasa_adjust[2], s.d_adjust = a_sasa_adjust[3])

# Cytosine p-values
c_sasa_df <- sasa_df %>% filter(don_resn == "C")
c_sasa_pvals <- data.frame(don_resn = "C",
                           n.s = wilcox.test(pull(filter(c_sasa_df, type == "Non"), SASA),
                                             pull(filter(c_sasa_df, type == "Single"), SASA))$p.value,
                           n.d = wilcox.test(pull(filter(c_sasa_df, type == "Non"), SASA),
                                             pull(filter(c_sasa_df, type == "Dual"), SASA))$p.value,
                           s.d = wilcox.test(pull(filter(c_sasa_df, type == "Single"), SASA),
                                             pull(filter(c_sasa_df, type == "Dual"), SASA))$p.value)
c_sasa_adjust <- p.adjust(as.numeric(c_sasa_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
c_sasa_pvals <-
  cbind(c_sasa_pvals, n.s_adjust = c_sasa_adjust[1], n.d_adjust = c_sasa_adjust[2], s.d_adjust = c_sasa_adjust[3])

# Guanine p-values
g_sasa_df <- sasa_df %>% filter(don_resn == "G")
g_sasa_pvals <- data.frame(don_resn = "G",
                           n.s = wilcox.test(pull(filter(g_sasa_df, type == "Non"), SASA),
                                             pull(filter(g_sasa_df, type == "Single"), SASA))$p.value,
                           n.d = wilcox.test(pull(filter(g_sasa_df, type == "Non"), SASA),
                                             pull(filter(g_sasa_df, type == "Dual"), SASA))$p.value,
                           s.d = wilcox.test(pull(filter(g_sasa_df, type == "Single"), SASA),
                                             pull(filter(g_sasa_df, type == "Dual"), SASA))$p.value)
g_sasa_adjust <- p.adjust(as.numeric(g_sasa_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
g_sasa_pvals <-
  cbind(g_sasa_pvals, n.s_adjust = g_sasa_adjust[1], n.d_adjust = g_sasa_adjust[2], s.d_adjust = g_sasa_adjust[3])

# Write SASA p-values to a csv file.
sasa_pvals <- rbind(a_sasa_pvals, c_sasa_pvals, g_sasa_pvals)
write.csv(sasa_pvals, file = snakemake@output[["sasa_pvals"]], row.names = FALSE)

#### HEAVY ATOM DENSITY PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
density_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
density_df[density_df$type == 0, "type"] <- "Non"
density_df[density_df$type == 1, "type"] <- "Single"
density_df[density_df$type == 2, "type"] <- "Dual"
density_df$type <- factor(density_df$type, levels = c("Non", "Single", "Dual"))

# Create a column that contains the donor atom name along with the donor residue name.
density_df <- density_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Add density columns.
volume_1 <- (4/3)*pi*snakemake@config[["count_dist_1"]]**3
volume_2 <- (4/3)*pi*snakemake@config[["count_dist_2"]]**3
volume_diff <- volume_2 - volume_1
density_df$density_1 <- density_df$count_1/volume_1
density_df$density_2 <- (density_df$count_2-density_df$count_1)/volume_diff

# Pivot the density_1 and density_2 columns to place the values in a single column.
density_pivot_df <- density_df %>%
  pivot_longer(cols = c(density_1, density_2), names_to = "density_type", values_to = "density_values")

# Change the density_type naming.
density_pivot_df[density_pivot_df$density_type == "density_1", "density_type"] <- "ROI 1"
density_pivot_df[density_pivot_df$density_type == "density_2", "density_type"] <- "ROI 2"

# Create the box plots with all data (outliers included).
density_plot <- density_pivot_df %>% ggplot(aes(x=type, y=density_values*1000, fill=density_type)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.stroke = 0, outlier.size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c(color_palette[3], color_palette[6])) +
  xlab("Type of Donating Amine") +
  ylab("Density (heavy-atoms/nm\ub3)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid.major.x = element_blank()) +
  facet_nested_wrap(vars(density_type, don_label), nrow = 1, scales = "fixed")

# Create the box plots without outliers.
density_no_outiers_plot <- density_pivot_df %>% ggplot(aes(x=type, y=density_values*1000, fill=density_type)) +
  geom_boxplot(outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = c(color_palette[3], color_palette[6])) +
  xlab("Type of Donating Amine") +
  ylab("Density (heavy-atoms/nm\ub3)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid.major.x = element_blank()) +
  facet_nested_wrap(vars(density_type, don_label), nrow = 1, scales = "fixed")

# Write the plots.
ggsave(snakemake@output[["density"]], plot = density_plot, width = 6, height = 3, units = "in", scale = 1)
ggsave(snakemake@output[["density_no_outiers"]],
       plot = density_no_outiers_plot, width = 6, height = 3, units = "in", scale = 1)

# Calculate quantile values for both densities.
density_stats <- density_pivot_df %>% summarise(n = n(),
                                                med = median(density_values),
                                                iqr = IQR(density_values),
                                                q1 = quantile(density_values, 0.25),
                                                q3 = quantile(density_values, 0.75),
                                                .by = c(type, don_label, density_type))

# Write density statistics to a csv file.
write.csv(density_stats, file = snakemake@output[["density_stats"]], row.names = FALSE)

#### HEAVY ATOM DENSITY WILCOXON RANK SUM TESTS ####

# Adenine ROI 1 p-values
a_roi_1_df <- density_pivot_df %>% filter(don_resn == "A", density_type == "ROI 1")
a_roi_1_pvals <- data.frame(don_resn = "A", density_type = "ROI 1",
                            n.s = wilcox.test(pull(filter(a_roi_1_df, type == "Non"), density_values),
                                              pull(filter(a_roi_1_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(a_roi_1_df, type == "Non"), density_values),
                                              pull(filter(a_roi_1_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(a_roi_1_df, type == "Single"), density_values),
                                              pull(filter(a_roi_1_df, type == "Dual"), density_values))$p.value)
a_roi_1_adjust <- p.adjust(as.numeric(a_roi_1_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
a_roi_1_pvals <-
  cbind(a_roi_1_pvals, n.s_adjust = a_roi_1_adjust[1], n.d_adjust = a_roi_1_adjust[2], s.d_adjust = a_roi_1_adjust[3])

# Cytosine ROI 1 p-values
c_roi_1_df <- density_pivot_df %>% filter(don_resn == "C", density_type == "ROI 1")
c_roi_1_pvals <- data.frame(don_resn = "C", density_type = "ROI 1",
                            n.s = wilcox.test(pull(filter(c_roi_1_df, type == "Non"), density_values),
                                              pull(filter(c_roi_1_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(c_roi_1_df, type == "Non"), density_values),
                                              pull(filter(c_roi_1_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(c_roi_1_df, type == "Single"), density_values),
                                              pull(filter(c_roi_1_df, type == "Dual"), density_values))$p.value)
c_roi_1_adjust <- p.adjust(as.numeric(c_roi_1_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
c_roi_1_pvals <-
  cbind(c_roi_1_pvals, n.s_adjust = c_roi_1_adjust[1], n.d_adjust = c_roi_1_adjust[2], s.d_adjust = c_roi_1_adjust[3])

# Guanine ROI 1 p-values
g_roi_1_df <- density_pivot_df %>% filter(don_resn == "G", density_type == "ROI 1")
g_roi_1_pvals <- data.frame(don_resn = "G", density_type = "ROI 1",
                            n.s = wilcox.test(pull(filter(g_roi_1_df, type == "Non"), density_values),
                                              pull(filter(g_roi_1_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(g_roi_1_df, type == "Non"), density_values),
                                              pull(filter(g_roi_1_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(g_roi_1_df, type == "Single"), density_values),
                                              pull(filter(g_roi_1_df, type == "Dual"), density_values))$p.value)
g_roi_1_adjust <- p.adjust(as.numeric(g_roi_1_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
g_roi_1_pvals <-
  cbind(g_roi_1_pvals, n.s_adjust = g_roi_1_adjust[1], n.d_adjust = g_roi_1_adjust[2], s.d_adjust = g_roi_1_adjust[3])

# Adenine ROI 2 p-values
a_roi_2_df <- density_pivot_df %>% filter(don_resn == "A", density_type == "ROI 2")
a_roi_2_pvals <- data.frame(don_resn = "A", density_type = "ROI 2",
                            n.s = wilcox.test(pull(filter(a_roi_2_df, type == "Non"), density_values),
                                              pull(filter(a_roi_2_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(a_roi_2_df, type == "Non"), density_values),
                                              pull(filter(a_roi_2_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(a_roi_2_df, type == "Single"), density_values),
                                              pull(filter(a_roi_2_df, type == "Dual"), density_values))$p.value)
a_roi_2_adjust <- p.adjust(as.numeric(a_roi_2_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
a_roi_2_pvals <-
  cbind(a_roi_2_pvals, n.s_adjust = a_roi_2_adjust[1], n.d_adjust = a_roi_2_adjust[2], s.d_adjust = a_roi_2_adjust[3])

# Cytosine ROI 2 p-values
c_roi_2_df <- density_pivot_df %>% filter(don_resn == "C", density_type == "ROI 2")
c_roi_2_pvals <- data.frame(don_resn = "C", density_type = "ROI 2",
                            n.s = wilcox.test(pull(filter(c_roi_2_df, type == "Non"), density_values),
                                              pull(filter(c_roi_2_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(c_roi_2_df, type == "Non"), density_values),
                                              pull(filter(c_roi_2_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(c_roi_2_df, type == "Single"), density_values),
                                              pull(filter(c_roi_2_df, type == "Dual"), density_values))$p.value)
c_roi_2_adjust <- p.adjust(as.numeric(c_roi_2_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
c_roi_2_pvals <-
  cbind(c_roi_2_pvals, n.s_adjust = c_roi_2_adjust[1], n.d_adjust = c_roi_2_adjust[2], s.d_adjust = c_roi_2_adjust[3])

# Guanine ROI 2 p-values
g_roi_2_df <- density_pivot_df %>% filter(don_resn == "G", density_type == "ROI 2")
g_roi_2_pvals <- data.frame(don_resn = "G", density_type = "ROI 2",
                            n.s = wilcox.test(pull(filter(g_roi_2_df, type == "Non"), density_values),
                                              pull(filter(g_roi_2_df, type == "Single"), density_values))$p.value,
                            n.d = wilcox.test(pull(filter(g_roi_2_df, type == "Non"), density_values),
                                              pull(filter(g_roi_2_df, type == "Dual"), density_values))$p.value,
                            s.d = wilcox.test(pull(filter(g_roi_2_df, type == "Single"), density_values),
                                              pull(filter(g_roi_2_df, type == "Dual"), density_values))$p.value)
g_roi_2_adjust <- p.adjust(as.numeric(g_roi_2_pvals[1, c("n.s", "n.d", "s.d")]), method = "BH")
g_roi_2_pvals <-
  cbind(g_roi_2_pvals, n.s_adjust = g_roi_2_adjust[1], n.d_adjust = g_roi_2_adjust[2], s.d_adjust = g_roi_2_adjust[3])

# Write density p-values to a csv file.
density_pvals <- rbind(a_roi_1_pvals, c_roi_1_pvals, g_roi_1_pvals, a_roi_2_pvals, c_roi_2_pvals, g_roi_2_pvals)
write.csv(density_pvals, file = snakemake@output[["density_pvals"]], row.names = FALSE)

#### PSEUDO-TORSION AND CORRESPONDING ACCEPTOR PAIR PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
pair_1_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G") & type == 2) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Translate the pseudo-torsions to range from 0 to 360 degrees.
pair_1_df[which(pair_1_df$eta >= 0), "eta_translated"] <- pair_1_df[which(pair_1_df$eta >= 0), "eta"]
pair_1_df[which(pair_1_df$eta < 0), "eta_translated"] <- pair_1_df[which(pair_1_df$eta < 0), "eta"] + 360
pair_1_df[which(pair_1_df$theta >= 0), "theta_translated"] <- pair_1_df[which(pair_1_df$theta >= 0), "theta"]
pair_1_df[which(pair_1_df$theta < 0), "theta_translated"] <- pair_1_df[which(pair_1_df$theta < 0), "theta"] + 360

# Filter for location 1.
pair_1_df <- pair_1_df %>% filter(eta_translated >= 43.2 & eta_translated < 72 &
                                    theta_translated >= 151.2 & theta_translated < 180)

# Filter for the top acceptor pair in location 1, only keep the necessary columns, and then write to a csv file.
pair_1_top_df <- pair_1_df %>% filter(don_resn == "A" & acc_pair_combined == "A(N7), N(NPO)")
pair_1_top_df <- pair_1_top_df[c("PDB", "model", "don_chain", "don_resi", "acc_pair_1_chain", "acc_pair_1_resi",
                                 "acc_pair_1_name", "acc_pair_2_chain", "acc_pair_2_resi", "acc_pair_2_name")]
write.csv(pair_1_top_df, snakemake@output[["location_1_residues"]], quote = FALSE, na = "NaN", row.names = FALSE)

# Create a column that contains the donor atom name along with the donor residue name.
pair_1_df <- pair_1_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Create a data frame with just acceptor pairs that belong to the same residue.
pair_1_same_df <- pair_1_df %>% filter(same_resi == "True")

# Calculate samples sizes and merge into the data frames.
pair_1_df <- merge(pair_1_df, summarise(pair_1_df, n_resn_pair = n(),
                                                  .by = c(don_label, acc_pair_combined)))
pair_1_same_df <- merge(pair_1_same_df, summarise(pair_1_same_df, n_resn_pair_same = n(),
                                                            .by = c(don_label, acc_pair_combined)))

# Only keep one row for each donor resn and acceptor pair combination.
pair_1_distinct_df <- pair_1_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)
pair_1_same_df <- pair_1_same_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)

# Add a column specifying the sample sizes for acceptor pairs that belong to the same residue.
pair_1_distinct_df <- merge(pair_1_distinct_df, pair_1_same_df[,c("don_label", "acc_pair_combined", "n_resn_pair_same")],
                        all.x = TRUE)

# Create a column that indicates the fraction of acceptor pairs that belong to the same residue. Fractions that would
# be 0 will be given a value of NA.
pair_1_distinct_df <- pair_1_distinct_df %>% mutate(fraction_label = NA)
pair_1_distinct_df[which(!is.na(pair_1_distinct_df$n_resn_pair_same)), "fraction_label"] <-
  paste(pair_1_distinct_df[which(!is.na(pair_1_distinct_df$n_resn_pair_same)), "n_resn_pair_same"], "/",
        pair_1_distinct_df[which(!is.na(pair_1_distinct_df$n_resn_pair_same)), "n_resn_pair"], sep = " ")

# Create a column that contains the combined acceptor pair divided by a forward slash instead of a comma.
pair_1_distinct_df <- pair_1_distinct_df %>%
  mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))

# Extract the rows from the combined data frame relevant for these plots.
pair_2_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G") & type == 2) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Translate the pseudo-torsions to range from 0 to 360 degrees.
pair_2_df[which(pair_2_df$eta >= 0), "eta_translated"] <- pair_2_df[which(pair_2_df$eta >= 0), "eta"]
pair_2_df[which(pair_2_df$eta < 0), "eta_translated"] <- pair_2_df[which(pair_2_df$eta < 0), "eta"] + 360
pair_2_df[which(pair_2_df$theta >= 0), "theta_translated"] <- pair_2_df[which(pair_2_df$theta >= 0), "theta"]
pair_2_df[which(pair_2_df$theta < 0), "theta_translated"] <- pair_2_df[which(pair_2_df$theta < 0), "theta"] + 360

# Filter for location 2.
pair_2_df <- pair_2_df %>% filter(eta_translated >= 295.2 & eta_translated < 324 &
                                    theta_translated >= 14.4 & theta_translated < 43.2)

# Filter for the top two acceptor pairs in location 2, only keep the necessary columns, and then write to a csv file.
pair_2_top_df <- pair_2_df %>% filter(don_resn == "G" & acc_pair_combined %in% c("N(NPO), U(O4)", "N(O5'), U(O4)"))
pair_2_top_df <- pair_2_top_df[c("PDB", "model", "don_chain", "don_resi", "acc_pair_1_chain", "acc_pair_1_resi",
                                 "acc_pair_1_name", "acc_pair_2_chain", "acc_pair_2_resi", "acc_pair_2_name")]
write.csv(pair_2_top_df, snakemake@output[["location_2_residues"]], quote = FALSE, na = "NaN", row.names = FALSE)

# Create a column that contains the donor atom name along with the donor residue name.
pair_2_df <- pair_2_df %>% mutate(don_label = paste(don_resn, "(", don_name, ")", sep = ""))

# Create a data frame with just acceptor pairs that belong to the same residue.
pair_2_same_df <- pair_2_df %>% filter(same_resi == "True")

# Calculate samples sizes and merge into the data frames.
pair_2_df <- merge(pair_2_df, summarise(pair_2_df, n_resn_pair = n(),
                                        .by = c(don_label, acc_pair_combined)))
pair_2_same_df <- merge(pair_2_same_df, summarise(pair_2_same_df, n_resn_pair_same = n(),
                                                  .by = c(don_label, acc_pair_combined)))

# Only keep one row for each donor resn and acceptor pair combination.
pair_2_distinct_df <- pair_2_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)
pair_2_same_df <- pair_2_same_df %>% distinct(don_label, acc_pair_combined, .keep_all = TRUE)

# Add a column specifying the sample sizes for acceptor pairs that belong to the same residue.
pair_2_distinct_df <- merge(pair_2_distinct_df, pair_2_same_df[,c("don_label", "acc_pair_combined", "n_resn_pair_same")],
                   all.x = TRUE)

# Create a column that indicates the fraction of acceptor pairs that belong to the same residue. Fractions that would
# be 0 will be given a value of NA.
pair_2_distinct_df <- pair_2_distinct_df %>% mutate(fraction_label = NA)
pair_2_distinct_df[which(!is.na(pair_2_distinct_df$n_resn_pair_same)), "fraction_label"] <-
  paste(pair_2_distinct_df[which(!is.na(pair_2_distinct_df$n_resn_pair_same)), "n_resn_pair_same"], "/",
        pair_2_distinct_df[which(!is.na(pair_2_distinct_df$n_resn_pair_same)), "n_resn_pair"], sep = " ")

# Create a column that contains the combined acceptor pair divided by a forward slash instead of a comma.
pair_2_distinct_df <- pair_2_distinct_df %>%
  mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))

# Combine the data frames.
pair_1_distinct_df["location"] <- "Location 1"
pair_2_distinct_df["location"] <- "Location 2"
combined_location_df <- rbind(pair_1_distinct_df, pair_2_distinct_df)

# Create a csv file with residues from Location 1 that are connected to the 5'-end of any residue from Location 2. This
# approach may not work if an insertion code is used for this or the downstream residue. If this is the case, print an
# error message, and create an empty csv file.
if (any(grepl("[a-zA-Z]", pair_1_df$don_resi)) | any(grepl("[a-zA-Z]", pair_2_df$don_resi))) {
  cat("Error: At least one residue within one of the pseudo-torsion locations contains an insertion code. The file
      named pt_neighbors.csv will be empty.")
  write.csv(data.frame(), "tables/pt_neighbors.csv", quote = FALSE, na = "NaN", row.names = FALSE)
} else {
  pair_1_df <- pair_1_df %>% mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))
  pair_2_df <- pair_2_df %>% mutate(acc_pair_combined_reformat = gsub(", ", "/", acc_pair_combined, fixed = TRUE))
  pair_1_df["location"] <- "Location 1"
  pair_2_df["location"] <- "Location 2"
  neighbors_df <- merge(pair_1_df %>% mutate(downstream_resi = as.character(as.numeric(don_resi) + 1)), pair_2_df,
                        by.x = c("downstream_resi", "don_chain", "eq_class_member"),
                        by.y = c("don_resi", "don_chain", "eq_class_member"), suffixes = c("",".y"))
  neighbors_df <- neighbors_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                                  "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                                  "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                                  "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                                  "PDB", "model", "eq_class_member",
                                  "eta", "theta", "eta_translated", "theta_translated", "chi")]
  write.csv(neighbors_df, snakemake@output[["pt_neighbors"]], quote = FALSE, na = "NaN", row.names = FALSE)
}

# Create csv files for the residues within Location 1 and Location 2.
write.csv(pair_1_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                       "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                       "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                       "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                       "PDB", "model", "eq_class_member",
                       "eta", "theta", "eta_translated", "theta_translated", "chi")],
          snakemake@output[["pt_location_1.csv"]], quote = FALSE, na = "NaN", row.names = FALSE)
write.csv(pair_2_df[,c("don_index", "don_name", "don_resn", "don_resi", "don_chain",
                       "acc_pair_1_name", "acc_pair_1_resi", "acc_pair_1_chain",
                       "acc_pair_2_name", "acc_pair_2_resi", "acc_pair_2_chain",
                       "don_label", "acc_pair_combined_reformat", "same_resi", "type",
                       "PDB", "model", "eq_class_member",
                       "eta", "theta", "eta_translated", "theta_translated", "chi")],
          snakemake@output[["pt_location_2.csv"]], quote = FALSE, na = "NaN", row.names = FALSE)

# Create the plots.
combined_location_plot <- combined_location_df %>%
  ggplot(aes(x=reorder_within(acc_pair_combined_reformat, n_resn_pair, list(don_label, location)),
             y=n_resn_pair, fill = don_label)) +
  geom_col(width = 0.8, color = "black") +
  geom_col_pattern(aes(y=n_resn_pair_same),
                   width = 0.8, color = "black", pattern = "stripe", pattern_fill = "black",
                   pattern_size = 0, pattern_angle = 45, pattern_spacing = 0.02, show.legend = FALSE) +
  geom_text(aes(y=n_resn_pair+max(n_resn_pair)*0.025, label=fraction_label),
            size = 8, size.unit = "pt", hjust = 0, vjust = 0, angle = 40) +
  xlab("Pair of Acceptor Atoms") +
  ylab("Count") +
  scale_fill_manual(values = c(color_palette[3], color_palette[8], color_palette[4]), name = "Donor") +
  scale_x_reordered(limits = function(x) rev(x)[1:4]) +
  scale_y_continuous(limits = c(0, 250)) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.key.size = unit(0.15, "in"), legend.key.spacing.y = unit(2, "pt")) +
  facet_wrap( ~ location, nrow = 1, scales = "free_x")

# The code for preparing the following plots is a bit more tedious because I wanted the plots to be the same size
# despite the legends having different widths due to the differing lengths of numbers in the legend keys.

# Extract the rows from the combined data frame relevant for these plots.
pt_df <- combined_df %>%
  filter(don_resn %in% c("A", "C", "G") & eta != -360 & theta != -360) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
pt_df[pt_df$type == 0, "type"] <- "Non"
pt_df[pt_df$type == 1, "type"] <- "Single"
pt_df[pt_df$type == 2, "type"] <- "Dual"
pt_df$type <- factor(pt_df$type, levels = c("Non", "Single", "Dual"))

# Translate the pseudo-torsions to range from 0 to 360 degrees.
pt_df[which(pt_df$eta >= 0), "eta_translated"] <- pt_df[which(pt_df$eta >= 0), "eta"]
pt_df[which(pt_df$eta < 0), "eta_translated"] <- pt_df[which(pt_df$eta < 0), "eta"] + 360
pt_df[which(pt_df$theta >= 0), "theta_translated"] <- pt_df[which(pt_df$theta >= 0), "theta"]
pt_df[which(pt_df$theta < 0), "theta_translated"] <- pt_df[which(pt_df$theta < 0), "theta"] + 360

# Set the width of the square bin.
bin_w <- 360/50

# Set the relative widths for the plots and legends.
grid_rel_w = c(5, 2)

# Set the boundaries of the helical region according to doi 10.1261/rna.079821.123.
eta_min <- 145
eta_max <- 190
theta_min <- 190
theta_max <- 245

# Create plots with the helical region and locations 1 and 2 outlined.
pt_no_plot_legend <- pt_df[which(pt_df$type == "Non"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(bin_w, bin_w)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.4, linetype=3, colour="red", fill=NA) + # helical region
  geom_rect(xmin = 43.2, xmax = 72, ymin = 151.2, ymax = 180,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 1
  geom_rect(xmin = 295.2, xmax = 324, ymin = 14.4, ymax = 43.2,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 2
  scale_fill_viridis(labels = comma, breaks = 10^seq(0, 3), name = "Count", transform = 'log10') +
  ggtitle("Non-Donating Amines") +
  xlab("\u03B7 (\ub0)") +
  ylab("\u03B8 (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(legend.key.width = unit(0.1, "in"), legend.key.height = unit(0.25, "in"),
        legend.margin = margin(t = 0, r = 0, b = 0.2, l = 0, "in"), plot.title = element_text(vjust = 1, hjust = 0.5),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0.1, "in"))
pt_single_plot_legend <- pt_df[which(pt_df$type == "Single"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(bin_w, bin_w)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.4, linetype=3, colour="red", fill=NA) + # helical region
  geom_rect(xmin = 43.2, xmax = 72, ymin = 151.2, ymax = 180,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 1
  geom_rect(xmin = 295.2, xmax = 324, ymin = 14.4, ymax = 43.2,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 2
  scale_fill_viridis(labels = comma, breaks = 10^seq(0, 4), name = "Count", transform = 'log10') +
  ggtitle("Single-Donating Amines") +
  xlab("\u03B7 (\ub0)") +
  ylab("\u03B8 (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(legend.key.width = unit(0.1, "in"), legend.key.height = unit(0.25, "in"),
        legend.margin = margin(t = 0, r = 0, b = 0.2, l = 0, "in"), plot.title = element_text(vjust = 1, hjust = 0.5),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0.1, "in"))
pt_dual_plot_legend <- pt_df[which(pt_df$type == "Dual"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(bin_w, bin_w)) +
  geom_rect(xmin = eta_min, xmax = eta_max, ymin = theta_min, ymax = theta_max,
            linewidth=0.4, linetype=3, colour="red", fill=NA) + # helical region
  geom_rect(xmin = 43.2, xmax = 72, ymin = 151.2, ymax = 180,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 1
  geom_rect(xmin = 295.2, xmax = 324, ymin = 14.4, ymax = 43.2,
            linewidth=0.4, linetype=1, colour="red", fill=NA) + # location 2
  scale_fill_viridis(labels = comma, breaks = 10^seq(0, 3), name = "Count", transform = 'log10') +
  ggtitle("Dual-Donating Amines") +
  xlab("\u03B7 (\ub0)") +
  ylab("\u03B8 (\ub0)") +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(legend.key.width = unit(0.1, "in"), legend.key.height = unit(0.25, "in"),
        legend.margin = margin(t = 0, r = 0, b = 0.2, l = 0, "in"), plot.title = element_text(vjust = 1, hjust = 0.5),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0.1, "in"))

# Extract the legends.
pt_no_legend <- get_legend(pt_no_plot_legend)
pt_single_legend <- get_legend(pt_single_plot_legend)
pt_dual_legend <- get_legend(pt_dual_plot_legend)

# Extract the plots.
pt_no_plot <- pt_no_plot_legend + theme(legend.position = "none")
pt_single_plot <- pt_single_plot_legend + theme(legend.position = "none")
pt_dual_plot <- pt_dual_plot_legend + theme(legend.position = "none")

# Combine the plots and legends.
pt_no_grid <- plot_grid(pt_no_plot, pt_no_legend, rel_widths = grid_rel_w)
pt_single_grid <- plot_grid(pt_single_plot, pt_single_legend, rel_widths = grid_rel_w)
pt_dual_grid <- plot_grid(pt_dual_plot, pt_dual_legend, rel_widths = grid_rel_w)

# Combine the plots.
pt_plots <- plot_grid(pt_no_grid, pt_single_grid, pt_dual_grid, combined_location_plot, nrow = 2) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Write the plots.
ggsave(snakemake@output[["pt"]], plot = pt_plots, width = 6.5, height = 5.5, units = "in")

# Write key data to a csv file.
pt_acc_pairs_df <-
  combined_location_df[c("location", "don_label", "acc_pair_combined_reformat", "n_resn_pair", "n_resn_pair_same")]
write.csv(pt_acc_pairs_df, file = snakemake@output[["pt_acc_pairs"]], row.names = FALSE)

#### CHI PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
chi_df <- combined_df %>% filter(don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
chi_df[chi_df$type == 0, "type"] <- "Non"
chi_df[chi_df$type == 1, "type"] <- "Single"
chi_df[chi_df$type == 2, "type"] <- "Dual"
chi_df$type <- factor(chi_df$type, levels = c("Non", "Single", "Dual"))

# Adjust the chi dihedrals such that values of -180 are instead 180.
chi_df <- chi_df %>% mutate(chi_adjusted = chi)
chi_df[which(chi_df$chi_adjusted == -180), "chi_adjusted"] <- 180

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

# Create the plots of the full 360 range for each don_resn and type combination.

# Find the maximum density for the full range.
max_full <- ceiling(max(chi_bins$density)*10^4)

# Adenine
chi_plot_full_1 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_2 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_3 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))

# Cytosine
chi_plot_full_4 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_5 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_6 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))

# Guanine
chi_plot_full_7 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_8 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))
chi_plot_full_9 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 0, ymax = max_full) +
  annotate("shadowtext", x = 75, y = 300,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 105, y = 300,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = "grey35") + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = pi) +
  scale_x_continuous(breaks = seq(-150, 180, 30), limits = c(-180, 180)) + scale_y_continuous(limits = c(0, max_full)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5, unit = "pt"))

# Convert the full plots into graphical objects.
grob_full_1 <- ggplotGrob(chi_plot_full_1)
grob_full_2 <- ggplotGrob(chi_plot_full_2)
grob_full_3 <- ggplotGrob(chi_plot_full_3)
grob_full_4 <- ggplotGrob(chi_plot_full_4)
grob_full_5 <- ggplotGrob(chi_plot_full_5)
grob_full_6 <- ggplotGrob(chi_plot_full_6)
grob_full_7 <- ggplotGrob(chi_plot_full_7)
grob_full_8 <- ggplotGrob(chi_plot_full_8)
grob_full_9 <- ggplotGrob(chi_plot_full_9)

# Combine the full plot grobs for each nucleobase type.
adenine_full <- cbind(grob_full_1, grob_full_2, grob_full_3)
cytosine_full <- cbind(grob_full_4, grob_full_5, grob_full_6)
guanine_full <- cbind(grob_full_7, grob_full_8, grob_full_9)

# Add x-axis, y-axis, and nucleobase type titles.
adenine_title <- ggdraw() + draw_label(expression(underline("Adenine")), size = 10, x = 0.53)
cytosine_title <- ggdraw() + draw_label(expression(underline("Cytosine")), size = 10, x = 0.53)
guanine_title <- ggdraw() + draw_label(expression(underline("Guanine")), size = 10, x = 0.53)
x_title <- ggdraw() + draw_label("\u03c7 (\ub0)", size = 10, x = 0.53)
y_title <- ggdraw() + draw_label(paste("Density (\ud7", "10\u207B\u2074)"), size = 10, angle = 90)

# Combine the rest of the grobs to create the completed set of plots.
full_no_y <- plot_grid(adenine_title, adenine_full, cytosine_title, cytosine_full, guanine_title, guanine_full,
                          x_title, ncol = 1, rel_heights = c(0.15, 1, 0.15, 1, 0.15, 1, 0.15))
full_plots <- plot_grid(y_title, full_no_y, nrow = 1, rel_widths = c(0.05, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Create a plot showing the combined counts over all amine donating and nucleobase types.

# Establish the cutoffs for the major and minor populations. They must be divisible by 10.
lower_cut <- -20
upper_cut <- 130

# Prepare a vector that contains the fill colors for the 36 bins in the resulting histograms.
chi_bin_colors <- c(rep("grey35", abs(-180-lower_cut)/10),
                    rep(color_palette[2], (upper_cut-lower_cut)/10),
                    rep("grey35", (180-upper_cut)/10))

# Create a plot showing all data combined.
chi_plot_combined <- chi_df %>% ggplot(aes(x = chi_adjusted)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -90, xmax = 90, ymin = 1, ymax = 10^5) +
  annotate("shadowtext", x = -120, y = 3*10^4,
           fontface = "bold.italic", label = 'anti', color = "white", bg.color = "black", bg.r = 0.075) +
  annotate("shadowtext", x = 0, y = 3*10^4,
           fontface = "bold.italic", label = 'syn', color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 135, y = 3*10^4,
           fontface = "bold.italic", label = 'anti', color = "white", bg.color = "black", bg.r = 0.075) +
  geom_histogram(aes(y=after_stat(count)), fill = chi_bin_colors, binwidth = 10, center = 5, show.legend = FALSE) +
  annotate("text", x = -105, y = 1*10^1, size = 4, label = 'major', fontface = "bold", color = "white") +
  annotate("text", x = 60, y = 1*10^1, size = 4, label = 'minor', fontface = "bold", color = "white") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 45)) +
  scale_y_continuous(transform = "log10", name = "Count") +
  xlab("\u03c7 (\ub0)") +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 0.4)

# Create the plots of a segment of the 360 range for each don_resn and type combination.

# Find the maximum density for the partial range.
max_partial <- ceiling(max(chi_bins[which(chi_bins$mids > -20 & chi_bins$mids < 130),]$density)*10^4)

# Adenine
chi_plot_partial_1 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_2 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_3 <- chi_bins[which(chi_bins$don_resn == "A" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))

# Cytosine
chi_plot_partial_4 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_5 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_6 <- chi_bins[which(chi_bins$don_resn == "C" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))

# Guanine
chi_plot_partial_7 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Non"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Non") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_8 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Single"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Single") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))
chi_plot_partial_9 <- chi_bins[which(chi_bins$don_resn == "G" & chi_bins$type == "Dual"),] %>%
  ggplot(aes(x=mids, y=density*10^4)) +
  annotate("rect", linetype = "blank", fill = color_palette[6], alpha = 0.3,
           xmin = -20, xmax = 90, ymin = 0, ymax = max_partial) +
  annotate("shadowtext", x = 0, y = 12,
           label = 'syn', size = 3, fontface = "bold.italic", color = color_palette[6], bg.color = "white") +
  annotate("shadowtext", x = 110, y = 10,
           label = 'anti', size = 3, fontface = "bold.italic", color = "white", bg.color = "black", bg.r = 0.075) +
  geom_col(fill = color_palette[2]) + ggtitle("Dual") +
  coord_radial(inner.radius = 0.3, expand = FALSE, start = -pi*(20/180), end = pi*(130/180)) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) + scale_y_continuous(limits = c(0, max_partial)) +
  xlab(element_blank()) + ylab(element_blank()) + theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5), panel.border = element_blank(),
        plot.margin = margin(t = 0, r = 5, b = 15, l = 5, unit = "pt"))

# Convert the partial plots into graphical objects.
grob_part_1 <- ggplotGrob(chi_plot_partial_1)
grob_part_2 <- ggplotGrob(chi_plot_partial_2)
grob_part_3 <- ggplotGrob(chi_plot_partial_3)
grob_part_4 <- ggplotGrob(chi_plot_partial_4)
grob_part_5 <- ggplotGrob(chi_plot_partial_5)
grob_part_6 <- ggplotGrob(chi_plot_partial_6)
grob_part_7 <- ggplotGrob(chi_plot_partial_7)
grob_part_8 <- ggplotGrob(chi_plot_partial_8)
grob_part_9 <- ggplotGrob(chi_plot_partial_9)

# Combine the partial plot grobs.
adenine_partial <- cbind(grob_part_1, grob_part_2, grob_part_3)
cytosine_partial <- cbind(grob_part_4, grob_part_5, grob_part_6)
guanine_partial <- cbind(grob_part_7, grob_part_8, grob_part_9)
partial_no_y <- plot_grid(adenine_title, adenine_partial, cytosine_title, cytosine_partial,
                          guanine_title, guanine_partial, x_title,
                          ncol = 1, rel_heights = c(0.15, 1, 0.15, 1, 0.15, 1, 0.15))
partial_plots <- plot_grid(y_title, partial_no_y, nrow = 1, rel_widths = c(0.1, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Create the plots.
ggsave(snakemake@output[["chi_full.png"]], plot = full_plots, width = 6.5, height = 7, units = "in", scale = 1)
ggsave(snakemake@output[["chi_combined.png"]], plot = chi_plot_combined,
  width = 3.75, height = 1.75, units = "in", scale = 1)
ggsave(snakemake@output[["chi_partial.png"]], plot = partial_plots, width = 4, height = 6.5, units = "in", scale = 1)

# Calculate the number of residues that belong to the major and minor conformations.
chi_df[!(chi_df["chi_adjusted"] > -20 & chi_df["chi_adjusted"] <= 130), "conformation"] <- "major"
chi_df[(chi_df["chi_adjusted"] > -20 & chi_df["chi_adjusted"] <= 130), "conformation"] <- "minor"
chi_stats <- chi_df %>% summarise(n = n(), .by = c(conformation, don_resn, type))

# Write the numbers to a csv file.
write.csv(chi_stats, file = snakemake@output[["chi_stats"]], row.names = FALSE)

#### CHI WATSON'S TWO-SAMPLE TESTS ####

# Create a new column based on chi_adjusted whose values are class circular.
chi_df <- chi_df %>% mutate(chi_circ = circular(chi_adjusted, units = "degrees", zero = -180))

# Start sinking the output into a text file.
sink(file = snakemake@output[["chi_pvals"]])

# Adenine p-values
a_chi_df <- chi_df %>% filter(don_resn == "A")
cat("A(N6) Non-Donating to Single-Donating Comparison:")
watson.two.test(pull(filter(a_chi_df, type == "Non"), chi_circ), pull(filter(a_chi_df, type == "Single"), chi_circ))
cat("A(N6) Non-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(a_chi_df, type == "Non"), chi_circ), pull(filter(a_chi_df, type == "Dual"), chi_circ))
cat("A(N6) Single-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(a_chi_df, type == "Single"), chi_circ), pull(filter(a_chi_df, type == "Dual"), chi_circ))

# Cytosine p-values
c_chi_df <- chi_df %>% filter(don_resn == "C")
cat("C(N4) Non-Donating to Single-Donating Comparison:")
watson.two.test(pull(filter(c_chi_df, type == "Non"), chi_circ), pull(filter(c_chi_df, type == "Single"), chi_circ))
cat("C(N4) Non-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(c_chi_df, type == "Non"), chi_circ), pull(filter(c_chi_df, type == "Dual"), chi_circ))
cat("C(N4) Single-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(c_chi_df, type == "Single"), chi_circ), pull(filter(c_chi_df, type == "Dual"), chi_circ))

# Guanine p-values
g_chi_df <- chi_df %>% filter(don_resn == "G")
cat("G(N2) Non-Donating to Single-Donating Comparison:")
watson.two.test(pull(filter(g_chi_df, type == "Non"), chi_circ), pull(filter(g_chi_df, type == "Single"), chi_circ))
cat("G(N2) Non-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(g_chi_df, type == "Non"), chi_circ), pull(filter(g_chi_df, type == "Dual"), chi_circ))
cat("G(N2) Single-Donating to Dual-Donating Comparison:")
watson.two.test(pull(filter(g_chi_df, type == "Single"), chi_circ), pull(filter(g_chi_df, type == "Dual"), chi_circ))

# Stop sinking the output into a text file.
sink(file = NULL)

# Remove the default Rplots.pdf generated.
file.remove('Rplots.pdf')
