# plot_pca.R
# Usage: Rscript 30_plot_pca.R input_matrix.tsv output_plot.pdf [highlight_sample] [psam_file] [color_by]

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
highlight_sample <- if(length(args) >= 3) args[3] else NA
psam_file <- if(length(args) >= 4) args[4] else NA
color_by <- if(length(args) >= 5) args[5] else "POP"  # Default to POP, can be "SUPERPOP"

# Load libraries
library(tidyverse)

# Read eigenvec file (plink2 format with #FID IID PC1 PC2 ...)
pca_df <- read_tsv(input_file, comment = "", show_col_types = FALSE) %>%
	rename(FID = 1, IID = 2) %>%  # Rename first two columns
	rename(Sample = IID)  # Use IID as Sample identifier

# If PSAM file provided, merge for POP and SUPERPOP info
if (!is.na(psam_file)) {
	# Read PSAM file (plink2 format with #FID IID PAT MAT SEX PHENO POP SUPERPOP...)
	psam_df <- read_tsv(psam_file, comment = "", show_col_types = FALSE) %>%
		rename(FID = 1, IID = 2) %>%
		rename(Sample = IID) %>%
		select(Sample, any_of(c("POP", "SUPERPOP")))
	
	pca_df <- pca_df %>%
		left_join(psam_df, by = "Sample")
}

# Add POP and SUPERPOP columns if missing
if (!"POP" %in% names(pca_df)) {
	pca_df <- pca_df %>%
		mutate(POP = NA_character_)
}

if (!"SUPERPOP" %in% names(pca_df)) {
	pca_df <- pca_df %>%
		mutate(SUPERPOP = NA_character_)
}

# Set the column to use for coloring
color_col <- color_by

# Set sample type (Main vs 1000G background)
pca_df <- pca_df %>%
	mutate(Type = if_else(Sample == highlight_sample & !is.na(highlight_sample), "Main", "Background"))

# Colour palette for Superpopulations (5 groups) - lighter shades
superpop_palette <- c(
	EUR = "#1976D2",  # Lighter blue
	AFR = "#FF6F00",  # Lighter orange
	EAS = "#388E3C",  # Lighter green
	SAS = "#D32F2F",  # Lighter red
	AMR = "#7B1FA2"   # Lighter purple
)

# Colour palette for Populations (26 populations from 1000G)
pop_palette <- c(
	# European
	CEU = "#1565C0", GBR = "#1976D2", FIN = "#1E88E5", IBS = "#2196F3", TSI = "#42A5F5",
	# African  
	YRI = "#E64A19", LWK = "#F4511E", GWD = "#FF5722", MSL = "#FF6F00", ESN = "#FF8A65", ACB = "#FFAB91", ASW = "#FFCCBC",
	# East Asian
	CHB = "#2E7D32", JPT = "#388E3C", CHS = "#43A047", CDX = "#4CAF50", KHV = "#66BB6A",
	# South Asian
	GIH = "#C62828", PJL = "#D32F2F", BEB = "#E53935", STU = "#F44336", ITU = "#EF5350",
	# American
	MXL = "#6A1B9A", PUR = "#7B1FA2", CLM = "#8E24AA", PEL = "#9C27B0"
)

main_colour <- "#FFD700"  # Gold for highlighted sample
background_alpha <- 0.4    # Semi-transparent for background
background_color <- "#999999"  # Medium grey for samples without metadata
sample_alpha <- 0.8        # Slightly transparent for our samples to show overlap
sample_size <- 2.8         # Size for our samples
background_size <- 2.0     # Size for background
highlight_size <- 3.2      # Larger for highlighted sample

# Shape assignments
background_shape <- 21  # Circle with border
sample_shape <- 24      # Triangle
highlight_shape <- 23   # Diamond

# Select the appropriate palette based on color_by
color_palette <- if(color_col == "SUPERPOP") superpop_palette else pop_palette

# Common theme for all plots
plot_theme <- theme_classic() +
	theme(
		legend.position = "right",
		plot.title = element_text(face = "bold", size = 16),
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12)
	)

# Function to create PCA plot
create_pca_plot <- function(df, pc_x, pc_y, title_base) {
	# Check if we have population data
	has_pop_data <- color_col %in% names(df) && any(!is.na(df[[color_col]]))
	
	# Identify our samples (those starting with T2T)
	df <- df %>%
		mutate(OurSample = str_detect(Sample, "^T2T"))
	
	# Create sample type for legend
	df <- df %>%
		mutate(SampleType = case_when(
			Type == "Main" ~ paste0("Highlighted: ", Sample),
			OurSample ~ "Study Samples",
			TRUE ~ "1000G Reference"
		))
	
	# Add highlighted sample to title only if a valid sample is specified
	title <- if (!is.na(highlight_sample) && highlight_sample != "" && highlight_sample != "NA") {
		paste0(title_base, " (Highlighted: ", highlight_sample, ")")
	} else {
		title_base
	}
	
	p <- ggplot(df, aes(x = .data[[pc_x]], y = .data[[pc_y]]))
	
	if (has_pop_data) {
		# Background samples with population colors (1000G reference)
		p <- p + geom_point(
			data = filter(df, Type == "Background" & !OurSample & !is.na(.data[[color_col]])),
			aes(color = .data[[color_col]], shape = SampleType), 
			alpha = background_alpha, 
			size = background_size,
			stroke = 0.5
		)
		
		# Our samples with population colors (more visible)
		p <- p + geom_point(
			data = filter(df, Type == "Background" & OurSample & !is.na(.data[[color_col]])),
			aes(color = .data[[color_col]], shape = SampleType), 
			alpha = sample_alpha, 
			size = sample_size,
			stroke = 0.8
		)
	}
	
	# Background samples without metadata
	if (any(df$Type == "Background" & !df$OurSample & (is.na(df[[color_col]]) | !has_pop_data))) {
		p <- p + geom_point(
			data = filter(df, Type == "Background" & !OurSample & (is.na(.data[[color_col]]) | !has_pop_data)),
			aes(shape = SampleType),
			color = background_color,
			alpha = background_alpha, 
			size = background_size,
			stroke = 0.5
		)
	}
	
	# Our samples without metadata
	if (any(df$Type == "Background" & df$OurSample & (is.na(df[[color_col]]) | !has_pop_data))) {
		p <- p + geom_point(
			data = filter(df, Type == "Background" & OurSample & (is.na(.data[[color_col]]) | !has_pop_data)),
			aes(shape = SampleType),
			color = "#555555",  # Medium grey
			alpha = sample_alpha, 
			size = sample_size,
			stroke = 0.8
		)
	}
	
	# Highlighted sample (if specified) - highest visibility
	if (any(df$Type == "Main")) {
		p <- p + geom_point(
			data = filter(df, Type == "Main"),
			aes(shape = SampleType),
			color = "grey20", 
			fill = main_colour,
			size = highlight_size,
			stroke = 0.8,
			alpha = 1.0
		)
	}
	
	# Add scales
	if (has_pop_data) {
		p <- p + scale_color_manual(
			name = color_col,
			values = color_palette, 
			na.value = background_color,
			drop = FALSE
		)
	}
	
	# Define shapes based on what sample types exist
	shape_values <- c(
		"1000G Reference" = background_shape,
		"Study Samples" = sample_shape
	)
	
	# Add highlighted sample shape if it exists
	if (any(df$Type == "Main")) {
		highlighted_label <- paste0("Highlighted: ", df$Sample[df$Type == "Main"][1])
		shape_values[highlighted_label] <- highlight_shape
	}
	
	p <- p + scale_shape_manual(
		name = "Sample Type",
		values = shape_values,
		breaks = names(shape_values)
	)
	
	p <- p + 
		labs(title = title, x = pc_x, y = pc_y) +
		plot_theme +
		guides(
			color = guide_legend(override.aes = list(size = 4, alpha = 1)),
			shape = guide_legend(override.aes = list(size = 4, alpha = 1))
		)
	
	return(p)
}

# Helper function to generate output filenames
get_output_name <- function(base_file, suffix) {
	# Remove any extension and add suffix
	base_name <- str_remove(base_file, "\\.[^.]+$")
	return(base_name)
}

base_name <- get_output_name(output_file, "")

# Main plot: PC1 vs PC2
gg_pc1_pc2 <- create_pca_plot(pca_df, "PC1", "PC2", "PCA Plot: PC1 vs PC2")
ggsave(paste0(base_name, ".pdf"), plot = gg_pc1_pc2, width = 7, height = 6)
ggsave(paste0(base_name, ".png"), plot = gg_pc1_pc2, width = 7, height = 6, dpi = 300)

# Additional plot: PC1 vs PC3
gg_pc1_pc3 <- create_pca_plot(pca_df, "PC1", "PC3", "PCA Plot: PC1 vs PC3")
ggsave(paste0(base_name, ".pc1_pc3.pdf"), plot = gg_pc1_pc3, width = 7, height = 6)
ggsave(paste0(base_name, ".pc1_pc3.png"), plot = gg_pc1_pc3, width = 7, height = 6, dpi = 300)

# Additional plot: PC2 vs PC3
gg_pc2_pc3 <- create_pca_plot(pca_df, "PC2", "PC3", "PCA Plot: PC2 vs PC3")
ggsave(paste0(base_name, ".pc2_pc3.pdf"), plot = gg_pc2_pc3, width = 7, height = 6)
ggsave(paste0(base_name, ".pc2_pc3.png"), plot = gg_pc2_pc3, width = 7, height = 6, dpi = 300)

# Grid views: All pairwise combinations of PC1-PC3 and PC1-PC4
library(patchwork)

# Create all pairwise plots for PC1-PC3 (without individual titles)
gg_grid_3_1_2 <- create_pca_plot(pca_df, "PC1", "PC2", "") + labs(title = NULL)
gg_grid_3_1_3 <- create_pca_plot(pca_df, "PC1", "PC3", "") + labs(title = NULL)
gg_grid_3_2_3 <- create_pca_plot(pca_df, "PC2", "PC3", "") + labs(title = NULL)

# Combine into grid layout for PC1-PC3 (1x3 grid)
gg_grid_pc3 <- (gg_grid_3_1_2 | gg_grid_3_1_3 | gg_grid_3_2_3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = if (!is.na(highlight_sample) && highlight_sample != "" && highlight_sample != "NA") {
      paste0("PCA Grid View: PC1-PC3 (Highlighted: ", highlight_sample, ")")
    } else {
      "PCA Grid View: PC1-PC3"
    },
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )

ggsave(paste0(base_name, ".grid_3pcs.pdf"), plot = gg_grid_pc3, width = 18, height = 6)
ggsave(paste0(base_name, ".grid_3pcs.png"), plot = gg_grid_pc3, width = 18, height = 6, dpi = 300)

# Create all pairwise plots for PC1-PC4 (without individual titles)
gg_grid_4_1_2 <- create_pca_plot(pca_df, "PC1", "PC2", "") + labs(title = NULL)
gg_grid_4_1_3 <- create_pca_plot(pca_df, "PC1", "PC3", "") + labs(title = NULL)
gg_grid_4_1_4 <- create_pca_plot(pca_df, "PC1", "PC4", "") + labs(title = NULL)
gg_grid_4_2_3 <- create_pca_plot(pca_df, "PC2", "PC3", "") + labs(title = NULL)
gg_grid_4_2_4 <- create_pca_plot(pca_df, "PC2", "PC4", "") + labs(title = NULL)
gg_grid_4_3_4 <- create_pca_plot(pca_df, "PC3", "PC4", "") + labs(title = NULL)

# Combine into grid layout for PC1-PC4 (3x2 grid)
gg_grid_pc4 <- (gg_grid_4_1_2 | gg_grid_4_1_3 | gg_grid_4_1_4) /
               (gg_grid_4_2_3 | gg_grid_4_2_4 | gg_grid_4_3_4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = if (!is.na(highlight_sample) && highlight_sample != "" && highlight_sample != "NA") {
      paste0("PCA Grid View: PC1-PC4 (Highlighted: ", highlight_sample, ")")
    } else {
      "PCA Grid View: PC1-PC4"
    },
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )

ggsave(paste0(base_name, ".grid_4pcs.pdf"), plot = gg_grid_pc4, width = 18, height = 12)
ggsave(paste0(base_name, ".grid_4pcs.png"), plot = gg_grid_pc4, width = 18, height = 12, dpi = 300)

# Variance explained plot (if eigenval file exists)
eigenval_file <- str_replace(input_file, "\\.eigenvec$", ".eigenval")
if (file.exists(eigenval_file)) {
	eigenvalues <- read_table(eigenval_file, col_names = "eigenval", show_col_types = FALSE)
	
	var_explained <- eigenvalues %>%
		mutate(
			PC = paste0("PC", row_number()),
			variance_pct = eigenval / sum(eigenval) * 100
		) %>%
		slice(1:10)
	
	gg_variance <- ggplot(var_explained, aes(x = fct_inorder(PC), y = variance_pct)) +
		geom_col(fill = "#1f77b4", alpha = 0.7) +
		labs(
			title = "Variance Explained by Principal Components", 
			x = "Principal Component", 
			y = "Variance Explained (%)"
		) +
		theme_classic() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
			plot.title = element_text(face = "bold", size = 16),
			axis.title = element_text(size = 14)
		)
	
	ggsave(paste0(base_name, ".variance_explained.pdf"), plot = gg_variance, width = 8, height = 5)
	ggsave(paste0(base_name, ".variance_explained.png"), plot = gg_variance, width = 8, height = 5, dpi = 300)
}
