#!/usr/bin/env Rscript

# Script to analyze RepeatMasker output files
# This script reads *.out files, summarizes repeat content,
# creates visualizations, and outputs BED files with repeat categories

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  #library(patchwork)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input RepeatMasker .out file", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="repeatmasker_analysis", 
              help="Output prefix for generated files [default= %default]", metavar="PREFIX"),
  make_option(c("-m", "--min_length"), type="integer", default=0,
              help="Minimum length of repeats to include [default= %default]"),
  make_option(c("-d", "--max_div"), type="numeric", default=100,
              help="Maximum divergence percentage to include [default= %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)) {
  stop("Input file must be specified with --input option")
}

# Read RepeatMasker .out file
message("Reading RepeatMasker output file: ", opt$input)
# Skip the header lines (first 3 lines)
repeatmasker_data <- read.table(opt$input, skip=3, stringsAsFactors=FALSE, fill=TRUE)

# Define column names based on RepeatMasker output format
colnames(repeatmasker_data) <- c(
  "score", "div", "del", "ins", "sequence", "begin", "end", "left", 
  "strand", "repeat_name", "class_family", "repeat_begin", "repeat_end", "repeat_left", "ID"
)

# Clean up the data
repeatmasker_data <- repeatmasker_data %>%
  filter(!is.na(score)) %>%
  mutate(
    length = end - begin + 1,
    class_family = gsub("^\\s+", "", class_family),
    repeat_class = ifelse(grepl("/", class_family), 
                         sub("/.*", "", class_family), 
                         class_family),
    repeat_family = ifelse(grepl("/", class_family), 
                          sub(".*/", "", class_family), 
                          "Unknown")
  )

# Create BED file for all repeats (unfiltered)
message("Creating unfiltered BED file...")
unfiltered_bed_data <- repeatmasker_data %>%
  select(sequence, begin, end, repeat_class, strand, repeat_family, repeat_name, div) %>%
  mutate(
    name = paste(repeat_class, repeat_family, sep=":"),
    rgb = case_when(
      repeat_class == "SINE" ~ "255,0,0",
      repeat_class == "LINE" ~ "0,255,0",
      repeat_class == "LTR" ~ "0,0,255",
      repeat_class == "DNA" ~ "255,255,0",
      repeat_class == "Simple_repeat" ~ "255,0,255",
      repeat_class == "Satellite" ~ "0,255,255",
      repeat_class == "Low_complexity" ~ "128,128,128",
      TRUE ~ "0,0,0"
    )
  )

# Write unfiltered BED file
write.table(
  unfiltered_bed_data %>% select(sequence, begin, end, name, repeat_name, strand, begin, end, rgb),
  paste0(opt$output, "_unfiltered.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Create a simplified unfiltered BED file with just major categories
unfiltered_simplified_bed <- repeatmasker_data %>%
  select(sequence, begin, end, repeat_class) %>%
  arrange(sequence, begin)

write.table(
  unfiltered_simplified_bed,
  paste0(opt$output, "_unfiltered_simplified.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Apply filters
filtered_data <- repeatmasker_data %>%
  filter(
    length >= opt$min_length,
    div <= opt$max_div
  )

message("Total repeats: ", nrow(repeatmasker_data))
message("Filtered repeats: ", nrow(filtered_data))

# Create summary statistics
# By class
class_summary <- filtered_data %>%
  group_by(repeat_class) %>%
  summarise(
    count = n(),
    total_length = sum(length),
    mean_length = mean(length),
    mean_div = mean(div)
  ) %>%
  arrange(desc(total_length))

# By family
family_summary <- filtered_data %>%
  group_by(repeat_class, repeat_family) %>%
  summarise(
    count = n(),
    total_length = sum(length),
    mean_length = mean(length),
    mean_div = mean(div)
  ) %>%
  arrange(repeat_class, desc(total_length))

# By sequence (chromosome)
sequence_summary <- filtered_data %>%
  group_by(sequence) %>%
  summarise(
    count = n(),
    total_length = sum(length),
    percent_repeat = NA,  # Will need sequence lengths to calculate this
    LINE = sum(length[repeat_class == "LINE"]),
    SINE = sum(length[repeat_class == "SINE"]),
    LTR = sum(length[repeat_class == "LTR"]),
    DNA = sum(length[repeat_class == "DNA"]),
    Simple_repeat = sum(length[repeat_class == "Simple_repeat"]),
    Satellite = sum(length[repeat_class == "Satellite"]),
    Low_complexity = sum(length[repeat_class == "Low_complexity"])
  )

# Create visualizations
message("Creating plots...")

# 1. Distribution of repeat classes
p1 <- ggplot(class_summary, aes(x = reorder(repeat_class, -total_length), y = total_length/1e6)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Repeat Class", y = "Total Length (Mb)", title = "Distribution of Repeat Classes")

# 2. Top 20 repeat families by length
top_families <- family_summary %>%
  arrange(desc(total_length)) %>%
  head(15)

p2 <- ggplot(top_families, aes(x = reorder(paste(repeat_class, repeat_family, sep=":"), -total_length), 
                             y = total_length/1e6, fill = repeat_class)) +
  geom_bar(stat = "identity", color = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Repeat Family", y = "Total Length (Mb)", title = "Top 15 Repeat Families") +
  scale_fill_viridis(discrete = TRUE)

# 3. Divergence distribution by repeat class
p3 <- ggplot(filtered_data, aes(x = div, fill = repeat_class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(x = "Divergence (%)", y = "Density", title = "Divergence Distribution by Repeat Class") +
  scale_fill_viridis(discrete = TRUE)

# 4. Length distribution by repeat class
p4 <- ggplot(filtered_data, aes(x = repeat_class, y = log10(length), fill = repeat_class)) +
  geom_boxplot(color = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Repeat Class", y = "Log10(Length)", title = "Length Distribution by Repeat Class") +
  scale_fill_viridis(discrete = TRUE)

# 5. Length histogram (non-log) by repeat class
p5 <- ggplot(filtered_data, aes(x = length, fill = repeat_class)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  facet_wrap(~repeat_class, scales = "free_y") +
  theme_classic() +
  labs(x = "Length (bp)", y = "Count", title = "Length Distribution by Repeat Class") +
  scale_fill_viridis(discrete = TRUE) +
  theme(legend.position = "none") +
  scale_x_log10()


# Save plots individually
ggsave(paste0(opt$output, "_class_summary.pdf"), p1, width = 14, height = 10)
ggsave(paste0(opt$output, "_family_summary.pdf"), p2, width = 14, height = 10)
ggsave(paste0(opt$output, "_divergence_distribution.pdf"), p3, width = 14, height = 10)
ggsave(paste0(opt$output, "_length_distribution.pdf"), p4, width = 14, height = 10)

ggsave(paste0(opt$output, "_class_summary.png"), p1, width = 14, height = 10, dpi = 300)
ggsave(paste0(opt$output, "_family_summary.png"), p2, width = 14, height = 10, dpi = 200)
ggsave(paste0(opt$output, "_divergence_distribution.png"), p3, width = 14, height = 10, dpi = 300)
ggsave(paste0(opt$output, "_length_distribution.png"), p4, width = 14, height = 10, dpi = 300)

# Save new plots
ggsave(paste0(opt$output, "_length_histogram.pdf"), p5, width = 14, height = 10)
ggsave(paste0(opt$output, "_length_histogram.png"), p5, width = 14, height = 10, dpi = 300)

# Write summary tables
write.csv(class_summary, paste0(opt$output, "_class_summary.csv"), row.names = FALSE)
write.csv(family_summary, paste0(opt$output, "_family_summary.csv"), row.names = FALSE)
write.csv(sequence_summary, paste0(opt$output, "_sequence_summary.csv"), row.names = FALSE)

# Create BED file with categories (filtered)
message("Creating filtered BED file...")
bed_data <- filtered_data %>%
  select(sequence, begin, end, repeat_class, strand, repeat_family, repeat_name, div) %>%
  mutate(
    name = paste(repeat_class, repeat_family, sep=":"),
    rgb = case_when(
      repeat_class == "SINE" ~ "255,0,0",
      repeat_class == "LINE" ~ "0,255,0",
      repeat_class == "LTR" ~ "0,0,255",
      repeat_class == "DNA" ~ "255,255,0",
      repeat_class == "Simple_repeat" ~ "255,0,255",
      repeat_class == "Satellite" ~ "0,255,255",
      repeat_class == "Low_complexity" ~ "128,128,128",
      TRUE ~ "0,0,0"
    )
  )

# Write filtered BED file
write.table(
  bed_data %>% select(sequence, begin, end, name, repeat_name, strand, begin, end, rgb),
  paste0(opt$output, "_filtered.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Create a simplified filtered BED file with just major categories
simplified_bed <- filtered_data %>%
  select(sequence, begin, end, repeat_class) %>%
  arrange(sequence, begin)

write.table(
  simplified_bed,
  paste0(opt$output, "_filtered_simplified.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

# Create a simplified bedfile with only Satellites
satellite_bed <- filtered_data %>%
  filter(repeat_class == "Satellite") %>%
  select(sequence, begin, end, repeat_family, repeat_name) %>%
  arrange(sequence, begin)

write.table(
  satellite_bed,
  paste0(opt$output, "_filtered_satellites.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)
# Create a bedfile for nucplot
message("Creating nucplot BED file...")
nucplot_bed <- filtered_data %>%
  select(sequence, begin, end, repeat_class) %>%
  mutate(plot = "plot") %>%
  arrange(sequence, begin)

write.table(
  nucplot_bed,
  paste0(opt$output, "_nucplot.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

message("Analysis complete. Output files saved with prefix: ", opt$output) 