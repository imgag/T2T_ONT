#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)

read_contig_info <- function(contig_file) {
  read_lines(contig_file) %>%
    str_subset("^#") %>%
    str_remove("^#\\s*") %>%
    read_tsv(col_names = c("contig", "length", "chromosome"))
}

read_reference_lengths <- function(fai_file) {
  read_lines(fai_file) %>%
    str_subset("^#") %>%
    str_remove("^#\\s*") %>%
    read_tsv(col_names = c("chromosome", "length", "offset", "linebases", "linewidth"))
}

read_bed_features <- function(bed_file) {
  read_lines(bed_file) %>%
    str_subset("^#") %>%
    str_remove("^#\\s*") %>%
    read_tsv(col_names = c("chromosome", "start", "end", "feature_type"))
}

create_karyogram <- function(contigs, ref_lengths, features, feature_colors) {
  contigs <- contigs %>%
    group_by(chromosome) %>%
    mutate(y_pos = row_number()) %>%
    ungroup()
  
  ggplot() +
    geom_rect(data = contigs,
              aes(xmin = 0, xmax = length,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4,
                  fill = contig),
              color = "black") +
    geom_vline(data = ref_lengths,
               aes(xintercept = length),
               linetype = "dashed", color = "red") +
    geom_rect(data = do.call(rbind, features),
              aes(xmin = start, xmax = end,
                  ymin = y_pos - 0.4, ymax = y_pos + 0.4,
                  fill = feature_type),
              alpha = 0.5) +
    facet_wrap(~chromosome, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = feature_colors) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size = 12),
          legend.position = "bottom") +
    labs(x = "Position (bp)",
         fill = "Feature")
}

main <- function(contig_file, fai_file, feature_files, output_file) {
  contigs <- read_contig_info(contig_file)
  ref_lengths <- read_reference_lengths(fai_file)
  features <- lapply(feature_files, read_bed_features)
  
  feature_colors <- c(
    "N_region" = "gray",
    "telomere" = "blue",
    "centromere" = "red"
  )
  
  p <- create_karyogram(contigs, ref_lengths, features, feature_colors)
  ggsave(output_file, p, width = 12, height = 8)
} 