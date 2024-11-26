library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)
library(tidyr)
library(tidytext)
library(svglite)
library(ggh4x)

#### H-BONDING MEASUREMENT PLOTS ####

# Creates data frames from the combined data.
combined_df <- read.csv(snakemake@input[["combined"]], header = TRUE,
                        na.strings = "NaN", comment.char = "#")

# Extract the rows from the combined data frame relevant for these plots.
pairs_df <- combined_df %>%
  filter(geom == 1 & don_resn %in% c("A", "C", "G")) %>%
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

#### PSEUDOTORSION PLOTS ####

# Extract the rows from the combined data frame relevant for these plots.
pt_df <- combined_df %>%
  filter(don_resn %in% c("A", "C", "G") &
           eta != -360 & theta != -360) %>%
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert donor type from numeric to string.
pt_df[pt_df$type == 0, "type"] <- "No"
pt_df[pt_df$type == 1, "type"] <- "Single"
pt_df[pt_df$type == 2, "type"] <- "Dual"
pt_df$type <- factor(pt_df$type, levels = c("No", "Single", "Dual"))

# Translate the pseudotorsions to range from 0 to 360 degrees.
pt_df[which(pt_df$eta >= 0), "eta_translated"] <-
  pt_df[which(pt_df$eta >= 0), "eta"]
pt_df[which(pt_df$eta < 0), "eta_translated"] <-
  pt_df[which(pt_df$eta < 0), "eta"] + 360
pt_df[which(pt_df$theta >= 0), "theta_translated"] <-
  pt_df[which(pt_df$theta >= 0), "theta"]
pt_df[which(pt_df$theta < 0), "theta_translated"] <-
  pt_df[which(pt_df$theta < 0), "theta"] + 360

# Create the plot.
pt_no <- pt_df[which(pt_df$type == "No"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("No") +
  xlab(expression(paste("Eta (\ub0)"))) +
  ylab(expression(paste("Theta (\ub0)"))) +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_single <- pt_df[which(pt_df$type == "Single"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Single") +
  xlab(expression(paste("Eta (\ub0)"))) +
  ylab(expression(paste("Theta (\ub0)"))) +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_dual <- pt_df[which(pt_df$type == "Dual"),] %>%
  ggplot(aes(x = eta_translated, y = theta_translated)) +
  geom_bin_2d(binwidth = c(360/100, 360/100)) +
  scale_fill_viridis(name = "Count") +
  ggtitle("Dual") +
  xlab(expression(paste("Eta (\ub0)"))) +
  ylab(expression(paste("Theta (\ub0)"))) +
  coord_fixed(ratio = 1, xlim = c(0, 360), ylim = c(0, 360)) +
  scale_x_continuous(breaks = seq(0, 360, 90)) +
  scale_y_continuous(breaks = seq(0, 360, 90)) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.125, "in"),
        legend.key.height = unit(0.2, "in"))
pt_plot <- grid.arrange(pt_no, pt_single, pt_dual, nrow = 3)

# Write the plot.
ggsave(snakemake@output[["pseudotorsion"]], plot = pt_plot, width = 3.25,
       height = 6, units = "in", scale = 1)

#### HEAVY ATOM DENSITY PLOTS ####

# Specify custom colors.
custom_greens <- RColorBrewer::brewer.pal(4, "Greens")

# Extract the rows from the combined data frame relevant for these plots.
density_df <- combined_df %>%
  filter(don_resn %in% c("A", "C", "G")) %>%
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
  filter(don_resn %in% c("A", "C", "G")) %>%
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
  filter(h_bond == 1 & don_resn %in% c("A", "C", "G")) %>%
  distinct(don_index, acc_index, eq_class_member, .keep_all = TRUE)

# Convert column to factor.
acc_id_df$type <- factor(acc_id_df$type, levels = c(0, 1, 2))

# Create a column with values that specify both the acceptor resn and name.
acc_id_df <-acc_id_df %>% mutate(acc_resn_name = paste(acc_resn, "(", acc_name,
                                                       ")", sep = ""))

# Combine values in the acc_resn_name column that involve common backbone acceptor atoms for RNA.
rna <- c("A", "C", "G", "U")
acc_id_df[acc_id_df$acc_resn %in% rna & acc_id_df$acc_name %in% c("OP1", "OP2"), "acc_resn_name"] <- "N(NPO)"
acc_id_df[acc_id_df$acc_resn %in% rna & acc_id_df$acc_name == "O2'", "acc_resn_name"] <- "N(O2')"
acc_id_df[acc_id_df$acc_resn %in% rna & acc_id_df$acc_name == "O3'", "acc_resn_name"] <- "N(O3')"
acc_id_df[acc_id_df$acc_resn %in% rna & acc_id_df$acc_name == "O4'", "acc_resn_name"] <- "N(O4')"
acc_id_df[acc_id_df$acc_resn %in% rna & acc_id_df$acc_name == "O5'", "acc_resn_name"] <- "N(O5')"

# Combine values in the acc_resn_name column that involve common backbone acceptor atoms for DNA.
dna <- c("DA", "DC", "DG", "DT")
acc_id_df[acc_id_df$acc_resn %in% dna & acc_id_df$acc_name %in% c("OP1", "OP2"), "acc_resn_name"] <- "DN(NPO)"
acc_id_df[acc_id_df$acc_resn %in% dna & acc_id_df$acc_name == "O3'", "acc_resn_name"] <- "DN(O3')"
acc_id_df[acc_id_df$acc_resn %in% dna & acc_id_df$acc_name == "O4'", "acc_resn_name"] <- "DN(O4')"
acc_id_df[acc_id_df$acc_resn %in% dna & acc_id_df$acc_name == "O5'", "acc_resn_name"] <- "DN(O5')"

# Combine values in the acc_resn_name column that involve common backbone acceptor atoms for amino acids.
amino_acids <- c("ALA", "ASP", "ASN", "ARG", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                 "PRO", "SER", "THR", "TRP", "TYR", "VAL")
acc_id_df[acc_id_df$acc_resn %in% amino_acids & acc_id_df$acc_name == "O", "acc_resn_name"] <- "AA(O)"

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

### CHI PLOT ###

# Extract the rows from the combined data frame relevant for these plots.
chi_df <- combined_df %>%
  filter(don_resn %in% c("A", "C", "G")) %>%
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
chi_bins <- data.frame(don_resn = character(), type = character(),
                       counts = integer(), density = double(), mids = integer())
for (split_idx in seq_along(chi_split)) {
  group_bin <- hist(chi_split[[split_idx]]$chi_adjusted,
                    breaks = seq(-180, 180, 10), plot = FALSE)
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
chi_plot_partial <- chi_bins %>% ggplot(aes(x=mids, y=density, fill=type)) +
  geom_col(show.legend = FALSE) +
  coord_radial(inner.radius = 0.3, expand = FALSE,
               start = -pi/8, end = pi*3/4) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) +
  scale_y_continuous(limits = c(0, 0.0018)) +
  xlab(expression(paste("\u03c7 (\ub0)"))) +
  ylab(element_blank()) +
  scale_fill_manual(values = custom_greens[2:4], name = "Type") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.3, "in"),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_nested_wrap(vars(don_resn, type), nrow = 3, scales = "fixed")

# Create the plots of a segment of the 360 range with the y-axis displayed.
chi_plot_partial_y <- chi_bins %>% ggplot(aes(x=mids, y=density, fill=type)) +
  geom_col(show.legend = FALSE) +
  coord_radial(inner.radius = 0.3, expand = FALSE,
               start = -pi/8, end = pi*3/4) +
  scale_x_continuous(limits = c(-20, 130), breaks = seq(0, 120, 30)) +
  scale_y_continuous(limits = c(0, 0.0018)) +
  xlab(expression(paste("\u03c7 (\ub0)"))) +
  ylab(element_blank()) +
  scale_fill_manual(values = custom_greens[2:4], name = "Type") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(t = 0, r = 0.1, b = 0, l = 0.3, "in"),
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10)) +
  facet_nested_wrap(vars(don_resn, type), nrow = 3, scales = "fixed")

# Create a plot showing all data combined.
chi_plot_combined <- chi_df %>% ggplot(aes(x=chi_adjusted)) +
  geom_histogram(aes(y=after_stat(count)), binwidth = 10, center = 5,
                 show.legend = FALSE) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 1,
           xend = -20, yend = 10^5) +
  annotate("segment", color = "red", linetype = "dashed", x = 130, y = 1,
           xend = 130, yend = 10^5) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 1,
           xend = 130, yend = 1) +
  annotate("segment", color = "red", linetype = "dashed", x = -20, y = 10^5,
           xend = 130, yend = 10^5) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 180, 45)) +
  scale_y_continuous(transform = "log10", name = "Count") +
  xlab(expression(paste("\u03c7 (\ub0)"))) +
  theme_bw(base_size = 10) +
  theme(aspect.ratio = 0.25)

# Write the plots.
ggsave(snakemake@output[["chi_partial"]], plot = chi_plot_partial, width = 3,
       height = 9, units = "in", scale = 1)
ggsave(snakemake@output[["chi_partial_y"]], plot = chi_plot_partial_y, width = 3,
       height = 9, units = "in", scale = 1)
ggsave(snakemake@output[["chi_combined"]], plot = chi_plot_combined, width = 6.5,
       height = 9, units = "in", scale = 1)

# Remove the default Rplots.pdf generated.
file.remove('Rplots.pdf')
