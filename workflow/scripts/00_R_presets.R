# R presets for plotting
# Source this file at the beginning of your R scripts

library(ggplot2)

# Recommendation (Nature Guidelines)
# Font: Helvetica
# Multi Part Figure Labels: 8pt bold, a b c d ..
# Max text size: 7pt
# Min text size: 5pt

# Custom ggplot theme based on Nature guidelines
# Usage: theme_set(t2t_theme()) or customize with t2t_theme(base_size = 7)
t2t_theme <- function(base_size = 6, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) +
        theme(
            # Text elements (in pt)
            plot.title = element_text(size = 8, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 7, hjust = 0),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            strip.text = element_text(size = 7, face = "bold"),
            
            # Background and grid
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.spacing = unit(0.5, "lines"),
            
            # Plot margins - add padding around data points
            plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt"),
            
            # Legend
            legend.position = "right",
            legend.key.size = unit(0.8, "lines")
        )
}

# Set as default theme
theme_set(t2t_theme())

# Define a color palettes
col_source <- c("HPRC" = "#9ba39f", "Samples" = "#bb1133")
col_haplotype <- c("Hap 1" = "#4292c6", "Hap 2" = "#ef6548")
