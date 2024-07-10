library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(svglite)

# Create a data frame from the data in the combined.csv file.
df <- read.csv(snakemake@input[["combined"]], header = TRUE, na.strings = "NaN")

# Convert columns to factors.
df$DOI <- factor(df$DOI, levels = c(0, 1))
df$geom <- factor(df$geom, levels = c(0, 1))
df$h_bond <- factor(df$h_bond, levels = c(0, 1))
df$type <- factor(df$type, levels = c(0, 1, 2))

#### HEAVY ATOM DENSITY PLOTS ####

# Extract the rows from the main data frame relevant for these plots.
density_df <- df %>% filter(DOI == 1) %>% 
  distinct(don_index, eq_class_member, .keep_all = TRUE)

# Convert DNA resn names to the corresponding RNA resn names
density_df[density_df$don_resn == "DA", "don_resn"] <- "A"
density_df[density_df$don_resn == "DC", "don_resn"] <- "C"
density_df[density_df$don_resn == "DG", "don_resn"] <- "G"

# Add density columns.
volume_1 <- (4/3)*pi*snakemake@config[["count_dist_1"]]**3
volume_2 <- (4/3)*pi*snakemake@config[["count_dist_2"]]**3
volume_diff <- volume_2 - volume_1
density_df$density_1 <- density_df$count_1/volume_1
density_df$density_2 <- (density_df$count_2-density_df$count_1)/volume_diff

# Calculate samples sizes.
density_summary <- density_df %>% summarise(n = n(), .by = c(type, don_resn))

# Specify the variables for the plots.
custom_colors <- RColorBrewer::brewer.pal(4, "Greens")
ratio_1 <- (3/(max(density_df$density_1) - min(density_df$density_1)))*1.5
ratio_2 <- (3/(max(density_df$density_2) - min(density_df$density_2)))*1.5
ylim_1 <- c(-0.004, max(density_df$density_1))
ylim_2 <- c(-0.004, max(density_df$density_2))

# Create the plots.
plot_density_1 <- density_df %>% ggplot(aes(x=type, y=density_1, fill=type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = density_summary,
            aes(x=type, y=-0.0025,
                label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 10, size.unit = "pt", vjust = 1) +
  coord_fixed(ratio = ratio_1, xlim = c(1, 3), ylim = ylim_1) +
  ggtitle("Atoms \u2264 " + snakemake@config[["count_dist_1"]] + " \uc5") +
  xlab("# of Amine H-Bonds") +
  ylab(expression(paste("Density (Atoms/\uc5\ub3)"))) +
  scale_fill_manual(values = custom_colors[2:4]) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_resn, nrow = 1)
plot_density_2 <- density_df %>% ggplot(aes(x=type, y=density_2, fill=type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  geom_text(data = density_summary,
            aes(x=type, y=-0.0025,
                label = paste("n = ", prettyNum(n, big.mark = ","), sep = "")),
            size = 10, size.unit = "pt", vjust = 1) +
  coord_fixed(ratio = ratio_2, xlim = c(1, 3), ylim = ylim_2) +
  ggtitle(snakemake@config[["count_dist_1"]] + " \uc5 \u003c Atoms \u2264 " +
            snakemake@config[["count_dist_2"]] + " \uc5") +
  xlab("# of Amine H-Bonds") +
  ylab(expression(paste("Density (Atoms/\uc5\ub3)"))) +
  scale_fill_manual(values = custom_colors[2:4]) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap( ~ don_resn, nrow = 1)
plot_density <- grid.arrange(plot_density_1, plot_density_2, nrow = 2)

# Write the plots.
ggsave(snakemake@output[["density"]], plot = plot_density, width = 6.5,
       height = 9, units = "in", scale = 1)
