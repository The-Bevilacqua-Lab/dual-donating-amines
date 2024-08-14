library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tidytext)
library(svglite)

#### H-BONDING MEASUREMENT PLOTS ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE,
                        na.strings = "NaN")

# Extract the rows from the combined data frame relevant for these plots.
pairs_df <- combined_df %>%
  filter(geom == 1 & DOI == 1 & don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, acc_index, eq_class_member, .keep_all = TRUE)

# Find the maximum fill value corresponding to the general region below the
# h_dist_max cutoff specified in the config file.
h_bond_region <- pairs_df[(pairs_df$h_acc_distance <= snakemake@config[["h_dist_max"]]),] %>%
  ggplot(aes(x = h_angle, y = h_acc_distance)) +
  geom_bin_2d(binwidth = c(120/100, 2.0/100))
max_value <- max(ggplot_build(h_bond_region)$data[[1]]$value)

# Create the plot.
pairs <- pairs_df %>% ggplot(aes(x = h_angle, y = h_acc_distance)) +
  geom_bin_2d(binwidth = c(120/100, 2.0/100)) +
  geom_segment(x=140, y=1.0, xend=140, yend=2.5, linewidth=0.4, linetype=2,
               colour="Red") +
  geom_segment(x=140, y=2.5, xend=180, yend=2.5, linewidth=0.4, linetype=2,
               colour="Red") +
  geom_segment(x=180, y=1.0, xend=180, yend=2.5, linewidth=0.4, linetype=2,
               colour="Red") +
  geom_segment(x=140, y=1.0, xend=180, yend=1.0, linewidth=0.4, linetype=2,
               colour="Red") +
  scale_fill_viridis(limits = c(1, max_value), name = "Count") +
  xlab(expression(paste("Angle (\ub0)"))) +
  ylab(expression(paste("Distance (\uc5)"))) +
  coord_fixed(ratio = 120/2.0, xlim = c(60, 180), ylim = c(1.0, 3.0)) +
  scale_x_continuous(breaks = seq(60, 180, 30)) +
  scale_y_continuous(breaks = seq(1.0, 3.0, 0.4)) +
  theme_classic(base_size = 10) +
  theme(legend.key.width = unit(0.125, "in"))

# Write the plot.
ggsave(snakemake@output[["pairs"]], plot = pairs, width = 3.25,
       height = 3.25, units = "in", scale = 1)

#### HEAVY ATOM DENSITY PLOTS ####

# Specify custom colors.
custom_greens <- RColorBrewer::brewer.pal(4, "Greens")

# Extract the rows from the combined data frame relevant for these plots.
density_df <- combined_df %>%
  filter(DOI == 1 & don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
density_df[density_df$type == 0, "type"] <- "No"
density_df[density_df$type == 1, "type"] <- "Single"
density_df[density_df$type == 2, "type"] <- "Dual"
density_df$type <- factor(density_df$type, levels = c("No",
                                                      "Single",
                                                      "Dual"))

# Add density columns.
volume_1 <- (4/3)*pi*snakemake@config[["count_dist_1"]]**3
volume_2 <- (4/3)*pi*snakemake@config[["count_dist_2"]]**3
volume_diff <- volume_2 - volume_1
density_df$density_1 <- density_df$count_1/volume_1
density_df$density_2 <- (density_df$count_2-density_df$count_1)/volume_diff

# Calculate samples sizes.
density_summary <- density_df %>% summarise(n = n(), .by = c(type, don_resn))

# Specify the variables for the plots.
ylim_1 <- c(-0.01, 0.08)
ylim_2 <- c(-0.01, 0.07)
ratio_1 <- (3/diff(ylim_1))*1.5
ratio_2 <- (3/diff(ylim_2))*1.5

# Create the plots.
density_1_plot <- density_df %>% ggplot(aes(x=type, y=density_1, fill=type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = density_summary,
            aes(x=type, y=-0.005,
                label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 10, size.unit = "pt", vjust = 1, angle = 30) +
  coord_fixed(ratio = ratio_1, xlim = c(1, 3), ylim = ylim_1) +
  ggtitle(paste("Atoms \u2264 ", snakemake@config[["count_dist_1"]], " \uc5",
                sep = "")) +
  xlab("Type of Amine H-Bonding") +
  ylab(expression(paste("Density (Atoms/\uc5\ub3)"))) +
  scale_fill_manual(values = custom_greens[2:4]) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_resn, nrow = 1)
density_2_plot <- density_df %>% ggplot(aes(x=type, y=density_2, fill=type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = density_summary,
            aes(x=type, y=-0.005,
                label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 10, size.unit = "pt", vjust = 1, angle = 30) +
  coord_fixed(ratio = ratio_2, xlim = c(1, 3), ylim = ylim_2) +
  ggtitle(paste(snakemake@config[["count_dist_1"]],
                " \uc5 \u003c Atoms \u2264 ",
                snakemake@config[["count_dist_2"]], " \uc5", sep = "")) +
  xlab("Type of Amine H-Bonding") +
  ylab(expression(paste("Density (Atoms/\uc5\ub3)"))) +
  scale_fill_manual(values = custom_greens[2:4]) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_resn, nrow = 1)
density_plot <- grid.arrange(density_1_plot, density_2_plot, nrow = 2)

# Write the plots.
ggsave(snakemake@output[["density"]], plot = density_plot, width = 6.5,
       height = 9, units = "in", scale = 1)

#### NUCLEOBASE IDENTITY PLOT ####

# Extract the rows from the combined data frame relevant for these plots.
nuc_id_df <- combined_df %>%
  filter(DOI == 1 & don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
nuc_id_df[nuc_id_df$type == 0, "type"] <- "No"
nuc_id_df[nuc_id_df$type == 1, "type"] <- "Single"
nuc_id_df[nuc_id_df$type == 2, "type"] <- "Dual"
nuc_id_df$type <- factor(nuc_id_df$type, levels = c("No",
                                                    "Single",
                                                    "Dual"))

# Calculate samples sizes and merge into dataframe.
nuc_id_df <- merge(nuc_id_df, summarise(nuc_id_df, n_resn = n(),
                                        .by = c(don_resn)))
nuc_id_df <- merge(nuc_id_df, summarise(nuc_id_df, n_resn_type = n(),
                                        .by = c(don_resn, type)))

# Only keep one row for each donor residue and donor type combination.
nuc_id_df <- nuc_id_df %>% distinct(don_resn, type, .keep_all = TRUE)

# Add a column that specifies the percent occurrence for each category.
nuc_id_df <-nuc_id_df %>% mutate(occurance = n_resn_type/sum(n_resn_type)*100)

# Create the plot.
nuc_id_plot <- nuc_id_df %>% ggplot(aes(x=type, y=occurance, fill=type)) +
  geom_col(width = 0.8, color = "black", show.legend = FALSE) +
  geom_text(aes(x=type, y=occurance+max(occurance)*0.05,
                label=paste(round(occurance, digits = 1),
                            "%", sep = "")), inherit.aes = FALSE) +
  xlab("Type of Amine H-Bonding") +
  ylab(expression(paste("Occurrence (%)"))) +
  scale_fill_manual(values = custom_greens[2:4]) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1.5) +
  facet_wrap( ~ don_resn, nrow = 1)

# Write the plots.
ggsave(snakemake@output[["nuc_id"]], plot = nuc_id_plot, width = 6.5,
       height = 9, units = "in", scale = 1)

#### ACCEPTOR IDENTITY PLOT ####

# Extract the rows from the combined data frame relevant for these plots.
acc_id_df <- combined_df %>%
  filter(DOI == 1 & h_bond == 1 & don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert column to factor.
acc_id_df$type <- factor(acc_id_df$type, levels = c(0, 1, 2))

# Create a column with values that specify both the acceptor resn and name.
acc_id_df <-acc_id_df %>% mutate(acc_resn_name = paste(acc_resn, "(", acc_name,
                                                       ")", sep = ""))

# Combine values in the acc_resn_name column that involve common backbone atoms.
acc_id_df[acc_id_df$acc_resn_name == "A(OP1)", "acc_resn_name"] <- "N(OP1)"
acc_id_df[acc_id_df$acc_resn_name == "A(OP2)", "acc_resn_name"] <- "N(OP2)"
acc_id_df[acc_id_df$acc_resn_name == "C(OP1)", "acc_resn_name"] <- "N(OP1)"
acc_id_df[acc_id_df$acc_resn_name == "C(OP2)", "acc_resn_name"] <- "N(OP2)"
acc_id_df[acc_id_df$acc_resn_name == "G(OP1)", "acc_resn_name"] <- "N(OP1)"
acc_id_df[acc_id_df$acc_resn_name == "G(OP2)", "acc_resn_name"] <- "N(OP2)"
acc_id_df[acc_id_df$acc_resn_name == "U(OP1)", "acc_resn_name"] <- "N(OP1)"
acc_id_df[acc_id_df$acc_resn_name == "U(OP2)", "acc_resn_name"] <- "N(OP2)"
acc_id_df[acc_id_df$acc_resn_name == "A(O2')", "acc_resn_name"] <- "N(O2')"
acc_id_df[acc_id_df$acc_resn_name == "C(O2')", "acc_resn_name"] <- "N(O2')"
acc_id_df[acc_id_df$acc_resn_name == "G(O2')", "acc_resn_name"] <- "N(O2')"
acc_id_df[acc_id_df$acc_resn_name == "U(O2')", "acc_resn_name"] <- "N(O2')"

# Calculate samples sizes and merge into data frame.
acc_id_df <- merge(acc_id_df, summarise(acc_id_df, n_resn_type = n(),
                                        .by = c(don_resn, acc_resn_name, type)))

# Only keep one row for each acceptor atom and donor type combination.
acc_id_df <- acc_id_df %>% distinct(don_resn, acc_resn_name, type,
                                    .keep_all = TRUE)

# Create the plot.
acc_id_1_plot <- acc_id_df %>% filter(type == 1) %>%
  ggplot(aes(x=reorder_within(acc_resn_name, n_resn_type, don_resn),
             y=n_resn_type, fill=type)) +
  geom_col(width = 0.8, color = "black", show.legend = FALSE) +
  ggtitle("Single H-Bonding Amines") +
  xlab("Acceptor Atom") +
  ylab(expression(paste("Count"))) +
  scale_x_reordered(limits = function(x) rev(x)[1:10]) +
  scale_fill_manual(values = custom_greens[3]) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap( ~ don_resn, nrow = 1, scales = "free")
acc_id_2_plot <- acc_id_df %>% filter(type == 2) %>%
  ggplot(aes(x=reorder_within(acc_resn_name, n_resn_type, don_resn),
             y=n_resn_type, fill=type)) +
  geom_col(width = 0.8, color = "black", show.legend = FALSE) +
  ggtitle("Dual H-Bonding Amines") +
  xlab("Acceptor Atom") +
  ylab(expression(paste("Count"))) +
  scale_x_reordered(limits = function(x) rev(x)[1:10]) +
  scale_fill_manual(values = custom_greens[4]) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_wrap( ~ don_resn, nrow = 1, scales = "free")
acc_id_plot <- grid.arrange(acc_id_1_plot, acc_id_2_plot, nrow = 2)

# Write the plots.
ggsave(snakemake@output[["acc_id"]], plot = acc_id_plot, width = 6.5,
       height = 9, units = "in", scale = 1)

### DISTANCE PLOT ###

# Creates data frames from the distance data.
distance_df <- read.csv(snakemake@input[["distances"]], header = TRUE, na.strings = "NaN")

# Order the acceptors of interest.
distance_df$acc_name_AOI <- factor(distance_df$acc_name_AOI, levels = c("N1", "O2", "N3", "O6", "N7"))

# Specify custom colors.
custom_blues <- RColorBrewer::brewer.pal(4, "Blues")
custom_oranges <- RColorBrewer::brewer.pal(4, "Oranges")

# Only keep data where the nucleobase of the acceptor of interest has an average b-factor below the
# distance_b_factor_cutoff specified in the config file.
distance_df <- distance_df[which(distance_df$b_factor < snakemake@config[["distance_b_factor_cutoff"]]),]

# Only keep distances that are greater than twice the standard H-O bond length
# (Gilli G.; Gilli P. The Nature of the Hydrogen Bond. 2009) to remove
# artifacts.
distance_df <- distance_df[which(distance_df$don_acc_distance_AOI > 2*0.983),]

# Convert donor type from numeric to string.
distance_df[distance_df$type == 0, "type"] <- "No"
distance_df[distance_df$type == 1, "type"] <- "Single"
distance_df[distance_df$type == 2, "type"] <- "Dual"
distance_df$type <- factor(distance_df$type, levels = c("No",
                                                      "Single",
                                                      "Dual"))

# Calculate samples sizes of non-overlapping interactions.
sample_sizes <- distance_df[which(distance_df$overlap == 0),] %>%
  summarise(n = n(), .by = c(acc_name_AOI, don_resn, type))

# Create the plot of non-overlapping interactions.
dist_plot <- distance_df[which(distance_df$overlap == 0),] %>%
  ggplot(aes(x=acc_name_AOI, y=don_acc_distance_AOI, fill=type)) +
  geom_point(color = "darkgrey", position = position_jitterdodge(), size = 0, stroke = 0.3, show.legend = FALSE) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_text(data = sample_sizes,
            aes(x=acc_name_AOI, y=snakemake@config[["don_dist_max"]]+0.1, group=type,
                label = prettyNum(n, big.mark = ",")),
            size = 8, size.unit = "pt", vjust = 1, angle = 30, position = position_dodge(width = 0.9)) +
  xlab("Acceptor") +
  ylab(expression(paste("Distance (\uc5)"))) +
  scale_fill_manual(values = custom_blues[2:4]) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 1.5) +
  facet_wrap( ~ don_resn, nrow = 1, scales = "free_x")

# Remove the default Rplots.pdf generated.
file.remove('Rplots.pdf')

# Write the plots.
ggsave(snakemake@output[["dist"]], plot = dist_plot, width = 6.5,
       height = 4.5, units = "in", scale = 1)
