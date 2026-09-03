#!/usr/bin/env Rscript

# Plot chromosome QC plots (karyograms)
# 1) Chromosome lengths vs ref length
# 2) Highlight possible errors and Genomic repeats (Telo, Centro)

# Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
    stop(paste(
        "Usage: Rscript plot_karyograms.R <fai> <seqinfo> <colors> <gap_stats>",
        "<nucflag_status> <rm_satellites> <telo> <output_prefix>",
        "\n\nArguments:",
        "\n  fai:             Reference genome .fai index file",
        "\n  seqinfo:         Assembly seqinfo.txt file",
        "\n  colors:          Colors mapping file (contig to chromosome)",
        "\n  gap_stats:       Gap statistics BED file",
        "\n  nucflag_status:  Nucflag misassembly BED file",
        "\n  rm_satellites:   RepeatMasker satellites BED file",
        "\n  telo:            Telomeric repeat windows CSV file",
        "\n  output_prefix:   Output file prefix (without extension)"
    ))
}

# Input files
fai <- args[1]
seqinfo <- args[2]
colors <- args[3]
gap_stats <- args[4]
nucflag_status <- args[5]
rm <- args[6]
telo <- args[7]
output_prefix <- args[8]

# Create output directory if it doesn't exist
output_dir <- dirname(output_prefix)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Load and prepare data
# ============================================================================

# Load ref chromosome lengths from Fasta Index
dt_fai <- read_tsv(fai, 
    col_names = c("chromosome", "length", "offset", "linebases", "linewidth", "qualoffset"),
    col_types = "cnnnnn",
    show_col_types = FALSE) %>%
    filter(chromosome != "chrM")

# Order reference chromosomes
chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")

dt_fai <- dt_fai %>%
    mutate(chromosome = factor(chromosome, levels = chromosome_order)) %>%
    filter(!is.na(chromosome))

# Load colors
dt_colors <- read_tsv(colors, col_types = "ccc", show_col_types = FALSE)

# Load assembly lengths (seqinfo)
dt_l <- read_tsv(seqinfo, col_names = c("contig", "length", "n_N"), 
                 col_types = "cnn", show_col_types = FALSE)

# Merge assembly lengths and colours
dt <- inner_join(dt_l, dt_colors, by = "contig") %>% 
    separate(col = contig, into = c("hap", NA), sep = "-", remove = FALSE) %>%
    mutate(hap = str_replace_all(hap, "lotype", "")) %>%
    mutate(chromosome = factor(chromosome, levels = chromosome_order)) %>%
    group_by(chromosome, hap) %>%
    arrange(desc(length), .by_group = TRUE) %>%
    mutate(
        contig_start_on_hap_bar = cumsum(lag(length, default = 0))
    ) %>%
    ungroup() %>%
    arrange(chromosome, hap, desc(length))

hap_plot_order <- rev(sort(unique(dt$hap)))

# ============================================================================
# Helper function to join BED files with assembly data
# ============================================================================

join_bed <- function(bed, dt, hap_levels_for_plot = hap_plot_order, 
                     feature_height = 0.8, feature_y_offset = -0.1) {
    contig_info_for_join <- dt %>%
        select(contig, chromosome, hap, contig_start_on_hap_bar) %>%
        distinct(contig, .keep_all = TRUE)

    features_to_plot <- bed %>%
        inner_join(contig_info_for_join, by = "contig") %>%
        mutate(
            abs_feature_start = contig_start_on_hap_bar + start,
            abs_feature_end = contig_start_on_hap_bar + end,
            y_center = abs_feature_start + (abs_feature_end - abs_feature_start) / 2,
            y_spread = abs_feature_end - abs_feature_start
        )

    return(features_to_plot)
}

# ============================================================================
# Load feature BED files
# ============================================================================

# Load gaps
dt_gaps <- read_tsv(gap_stats, 
                    col_names = c("contig", "start", "end", "label", "length"), 
                    col_types = "cddcd", show_col_types = FALSE)
dt_gaps_plot <- join_bed(dt_gaps, dt)

# Load nucflag
dt_nucflag <- read_tsv(nucflag_status, 
                       col_names = c("contig", "start", "end", "label", "length"), 
                       col_types = "cddcd", show_col_types = FALSE)
dt_nucflag_plot <- join_bed(dt_nucflag, dt)

# Load satellites
dt_satellites <- read_tsv(rm, 
                          col_names = c("contig", "start", "end", "family", "name"), 
                          col_types = "cddcc", show_col_types = FALSE)
dt_satellites_plot <- join_bed(dt_satellites, dt)

# Load telomeres from telomeric repeat windows
dt_telo <- read_csv(telo, col_types = "cddd", show_col_types = FALSE) %>%
    separate(id, sep = "_", into = c("contig", "location")) %>%
    group_by(contig) %>%
    mutate(repeat_number = max(forward_repeat_number, reverse_repeat_number)) %>%
    ungroup() %>%
    filter(repeat_number > 10) %>%
    select(contig, location, repeat_number) %>%
    inner_join(dt_l, by = "contig") %>%
    mutate(start = ifelse(location == "Start", 0, length - repeat_number * 5)) %>%
    mutate(end = ifelse(location == "Start", repeat_number * 5, length)) %>%
    select(contig, start, end)

dt_telo_plot <- join_bed(dt_telo, dt)

# ============================================================================
# Generate plots
# ============================================================================

# Plot 1: With Gaps, all contigs including unassigned
p1 <- ggplot(dt, aes(x = hap, y = length)) +
    geom_hline(data = dt_fai,
               aes(yintercept = length), 
               lwd = 1, lty = 1, colour = "blue") +
    geom_bar(stat = "identity", col = "grey20", fill = "grey90", width = 0.7) + 
    geom_tile(data = dt_satellites_plot,
              aes(y = y_center, height = y_spread),
              fill = "grey20", linewidth = 1, alpha = 0.7, width = 0.7) +
    geom_tile(data = dt_telo_plot,
              aes(y = y_center, height = y_spread),
              fill = "darkgreen", colour = "darkgreen", linewidth = 1.5, alpha = 0.7, width = 0.9) +
    geom_tile(data = dt_gaps_plot,
              aes(y = y_center, height = y_spread),
              fill = "red", colour = "red", linewidth = 0.8, alpha = 0.9, width = 0.8) +
    coord_flip() +
    facet_grid(rows = vars(chromosome)) +
    theme_classic() +
    guides(fill = "none") +
    scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), 
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(limits = hap_plot_order) +
    expand_limits(y = -1e6) +
    labs(y = "Assembled Length (bp)", x = "") +
    theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.background = element_blank(),
        panel.spacing.y = unit(0.2, "lines"),
        legend.position = "none",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )

ggsave(plot = p1, filename = paste0(output_prefix, ".karyogram.gaps.with_unassigned.png"), 
       height = 15, width = 8)

# Plot 2: With Gaps, Clean (no unassigned, only chr scale contigs)
p2 <- ggplot(dt %>%
                 filter(length >= 10e6) %>%
                 filter(hap != "unassigned"),
             aes(x = hap, y = length)) +
    geom_hline(data = dt_fai,
               aes(yintercept = length), 
               lwd = 1, lty = 1, colour = "blue") +
    geom_bar(stat = "identity", col = "grey20", fill = "grey90", width = 0.7) + 
    geom_tile(data = dt_satellites_plot,
              aes(y = y_center, height = y_spread),
              fill = "grey20", linewidth = 1, alpha = 0.7, width = 0.7) +
    geom_tile(data = dt_telo_plot,
              aes(y = y_center, height = y_spread),
              fill = "darkgreen", colour = "darkgreen", linewidth = 1.5, alpha = 0.7, width = 0.9) +
    geom_tile(data = dt_gaps_plot,
              aes(y = y_center, height = y_spread),
              fill = "red", colour = "red", linewidth = 0.8, alpha = 0.9, width = 0.8) +
    coord_flip() +
    facet_grid(rows = vars(chromosome)) +
    theme_classic() +
    guides(fill = "none") +
    scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), 
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(limits = c("hap2", "hap1")) +
    expand_limits(y = -1e6) +
    labs(y = "Assembled Length (bp)", x = "") +
    theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.background = element_blank(),
        panel.spacing.y = unit(0.2, "lines"),
        legend.position = "none",
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )

ggsave(plot = p2, filename = paste0(output_prefix, ".karyogram.gaps.clean.png"), 
       height = 15, width = 8)

# Plot 3: With Nucflag misassemblies, Clean (no unassigned, only chr scale contigs)
p3 <- ggplot(dt %>%
                 filter(length >= 10e6) %>%
                 filter(hap != "unassigned"),
             aes(x = hap, y = length)) +
    geom_hline(data = dt_fai,
               aes(yintercept = length), 
               lwd = 1, lty = 1, colour = "blue") +
    geom_bar(stat = "identity", col = "grey20", fill = "grey90", width = 0.7) + 
    geom_tile(data = dt_satellites_plot,
              aes(y = y_center, height = y_spread),
              fill = "grey20", linewidth = 1, alpha = 0.7, width = 0.7) +
    geom_tile(data = dt_telo_plot,
              aes(y = y_center, height = y_spread),
              fill = "darkgreen", colour = "darkgreen", linewidth = 1.5, alpha = 0.7, width = 0.9) +
    geom_tile(data = dt_nucflag_plot %>% filter(label != "HET"),
              aes(y = y_center, height = y_spread, fill = label, colour = label),
              linewidth = 1, alpha = 0.7, width = 0.7) +
    coord_flip() +
    facet_grid(rows = vars(chromosome)) +
    theme_classic() +
    scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), 
                       expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(limits = c("hap2", "hap1")) +
    expand_limits(y = -1e6) +
    labs(y = "Assembled Length (bp)", x = "") +
    theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
        strip.background = element_blank(),
        panel.spacing.y = unit(0.2, "lines"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()
    )

ggsave(plot = p3, filename = paste0(output_prefix, ".karyogram.nucflag.clean.png"), 
       height = 15, width = 8)

# Optional: Generate summary plots for gap lengths and nucflag types
dt_gaps_length <- dt_gaps %>%
    inner_join(dt, by = "contig") %>%
    mutate(contig_type = ifelse(length.y > 10e6, "full chrom", "small"))

p0 <- ggplot(dt_gaps_length, aes(x = length.x, fill = contig_type)) +
    geom_histogram(colour = "grey20") +
    theme_classic() +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8)) +
    labs(x = "Length (bp)", y = "Gaps in assembly", fill = "Scaffold type")

ggsave(plot = p0, filename = paste0(output_prefix, ".gap_lengths.png"), 
       height = 3, width = 4)

dt_nucflag_l <- dt_nucflag %>%
    inner_join(dt, by = "contig") %>%
    mutate(contig_type = ifelse(length > 10e6, "small", "full chrom")) %>%
    mutate(contig_type = factor(contig_type, levels = c("small", "full chrom")))

p01 <- ggplot(dt_nucflag_l, aes(x = label, fill = contig_type)) +
    geom_histogram(stat = "count", colour = "grey20") +
    theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.8),
          axis.text.x = element_text(size = 5)) +
    labs(x = "", y = "", fill = "Scaffold type")

ggsave(plot = p01, filename = paste0(output_prefix, ".nucflag_types.png"), 
       height = 3, width = 4)

cat("Karyogram plots generated successfully!\n")
cat("Output files:\n")
cat("  -", paste0(output_prefix, ".karyogram.gaps.with_unassigned.png\n"))
cat("  -", paste0(output_prefix, ".karyogram.gaps.clean.png\n"))
cat("  -", paste0(output_prefix, ".karyogram.nucflag.clean.png\n"))
cat("  -", paste0(output_prefix, ".gap_lengths.png\n"))
cat("  -", paste0(output_prefix, ".nucflag_types.png\n"))
