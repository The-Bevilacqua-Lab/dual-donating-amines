library(ggplot2)
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

# Create plots depicting the heavy atom density around the donors of interest.

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

# Combine the density columns by pivoting.
density_df <- density_df %>% 
  pivot_longer(c("density_1", "density_2"), 
               names_to = "density_set", values_to = "density")
density_df[density_df$density_set == "density_1", "density_set"] <- 
  "Atoms \u2264 5 \uc5"
density_df[density_df$density_set == "density_2", "density_set"] <- 
  "5 \uc5 \u003c Atoms \u2264 10 \uc5"

# Calculate samples sizes.
density_summary <- density_df %>% summarise(n = n(), .by = c(type, density_set))

# Specify the variables for the plots.
custom_colors <- RColorBrewer::brewer.pal(4, "Greens")
x_labels <- c(
  paste("0\nn = ", prettyNum(density_summary$n[1], big.mark = ","), sep = ""), 
  paste("1\nn = ", prettyNum(density_summary$n[3], big.mark = ","), sep = ""),
  paste("2\nn = ", prettyNum(density_summary$n[5], big.mark = ","), sep = ""))
ratio <- (3/(max(density_df$density) - min(density_df$density)))*1.5 
ylim <- c(min(density_df$density), max(density_df$density))

# Create the plots.
plot_density <- density_df %>% ggplot(aes(x=type, y=density, fill=type)) +
  geom_violin(show.legend = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0, outlier.size = 1, show.legend = FALSE) +
  coord_fixed(ratio = ratio, xlim = c(1, 3), ylim = ylim) +
  xlab(element_blank()) +
  ylab(expression(paste("Density (Atoms/\uc5\ub3)"))) +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(values = custom_colors[2:4]) +
  theme_bw(base_size = 10) +
  facet_wrap( ~ density_set + don_resn, nrow = 2)

# Write the plots.
ggsave(snakemake@output[["density"]], plot = plot_density, width = 6.5,
       height = 9, units = "in", scale = 1)
