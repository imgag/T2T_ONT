#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

#' Read and process contig lengths and chromosome assignments
#' @param contig_file Path to file with contig lengths and chromosome assignments
#' @return Dataframe with processed contig information
read_contig_info <- function(contig_file) {
  # Expected format: contig_name, length, chromosome
  read_tsv(contig_file, col_names = c("contig", "length", "chromosome"))
}

#' Read reference chromosome lengths from fai file
#' @param fai_file Path to reference genome fai file
#' @return Dataframe with chromosome lengths
read_reference_lengths <- function(fai_file) {
  # Expected format: chromosome, length, offset, linebases, linewidth
  read_tsv(fai_file, col_names = c("chromosome", "length", "offset", "linebases", "linewidth"))
}

#' Read BED files with genomic features
#' @param bed_file Path to BED file
#' @return Dataframe with feature information
read_bed_features <- function(bed_file) {
  # Expected format: chromosome, start, end, feature_type
  read_tsv(bed_file, col_names = c("chromosome", "start", "end", "feature_type"))
}

#' Create karyogram plot
#' @param contigs Contig information dataframe
#' @param ref_lengths Reference chromosome lengths dataframe
#' @param features List of dataframes with genomic features
#' @param feature_colors Named vector of colors for features
#' @return ggplot object
create_karyogram <- function(contigs, ref_lengths, features, feature_colors) {
  # Calculate y positions for contigs
  contigs <- contigs %>%
    group_by(chromosome) %>%
    mutate(y_pos = row_number()) %>%
    ungroup()
  
  # Create base plot
  ggplot() +
    # Plot contigs as rectangles
    geom_rect(data = contigs,
              aes(xmin = 0, xmax = length,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4,
                  fill = contig),
              color = "black") +
    # Add reference lengths as vertical lines
    geom_vline(data = ref_lengths,
               aes(xintercept = length),
               linetype = "dashed", color = "red") +
    # Add genomic features
    geom_rect(data = do.call(rbind, features),
              aes(xmin = start, xmax = end,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4,
                  fill = feature_type),
              alpha = 0.5) +
    # Facet by chromosome
    facet_wrap(~chromosome, scales = "free_y", ncol = 1) +
    # Customize appearance
    scale_fill_manual(values = feature_colors) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size = 12),
          legend.position = "bottom") +
    labs(x = "Position (bp)",
         fill = "Feature")
}

#' Main function to create karyogram
#' @param contig_file Path to contig information file
#' @param fai_file Path to reference genome fai file
#' @param feature_files List of paths to BED files with features
#' @param output_file Path to save the plot
main <- function(contig_file, fai_file, feature_files, output_file) {
  # Read input data
  contigs <- read_contig_info(contig_file)
  ref_lengths <- read_reference_lengths(fai_file)
  
  # Read feature files
  features <- lapply(feature_files, read_bed_features)
  
  # Define colors for features
  feature_colors <- c(
    "N_region" = "gray",
    "telomere" = "blue",
    "centromere" = "red"
  )
  
  # Create and save plot
  p <- create_karyogram(contigs, ref_lengths, features, feature_colors)
  ggsave(output_file, p, width = 12, height = 8)
}

# Example usage:
# main(
#   contig_file = "contigs.tsv",
#   fai_file = "reference.fai",
#   feature_files = c("n_regions.bed", "telomeres.bed", "centromeres.bed"),
#   output_file = "karyogram.pdf"
# ) 