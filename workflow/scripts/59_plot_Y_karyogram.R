#!/usr/bin/env Rscript
# Plot chrY assembly length and error tracks for all samples.
# Produces two figures:
#   <prefix>.Y_karyogram.pdf/png           — T2T contig + scaffold samples, contigs >1 Mb
#   <prefix>.Y_karyogram_fragmented.pdf/png — fragmented samples, all contigs
#
# Usage:
#   Rscript 59_plot_Y_karyogram.R <samples_yml> <output_prefix>
#
# Inputs (resolved relative to project root):
#   assembly/qc/phased_verkko/{sample}/colors.tsv
#   assembly/qc/phased_verkko/{sample}/both.mapped_T2T.paf
#   assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed
#   assembly/qc/phased_verkko/{sample}/T2T_contigs.both_motif.csv
#   assembly/qc/phased_verkko/{sample}/T2T_scaffolds.both.txt
#   analysis_other/repeatmasker/{sample}/rm_summary/{sample}_filtered_satellites.bed

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript 59_plot_Y_karyogram.R <samples_yml> <output_prefix>")
}

samples_yml   <- args[1]
output_prefix <- args[2]

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

# ── Constants ─────────────────────────────────────────────────────────────────
REF_Y_LENGTH <- 62460029  # CHM13v2.0 chrY
QC_DIR       <- "assembly/qc/phased_verkko"
RM_DIR       <- "analysis_other/repeatmasker"
BAR_WIDTH    <- 0.7
MIN_GAP_VIS  <- 100000

# ── Load sample list ──────────────────────────────────────────────────────────
samples <- readLines(samples_yml) |>
    trimws() |>
    (\(x) x[nchar(x) > 0 & !startsWith(x, "#") & x != "T2T10"])()

# ── 1. chrY contig names from colors.tsv ─────────────────────────────────────
y_colors <- map_dfr(samples, function(s) {
    f <- file.path(QC_DIR, s, "colors.tsv")
    if (!file.exists(f)) return(NULL)
    read_tsv(f, col_types = "ccc", show_col_types = FALSE) |>
        filter(chromosome == "chrY") |>
        mutate(sample = s)
})

if (nrow(y_colors) == 0) stop("No chrY contigs found for any sample.")

# ── 2. Contig lengths from PAF (col2 = query length) ─────────────────────────
y_lengths <- map_dfr(unique(y_colors$sample), function(s) {
    paf <- file.path(QC_DIR, s, "both.mapped_T2T.paf")
    if (!file.exists(paf)) return(NULL)
    contigs <- y_colors |> filter(sample == s) |> pull(contig)
    read_tsv(paf, col_names = FALSE, col_types = cols(.default = "c"),
             show_col_types = FALSE) |>
        filter(X1 %in% contigs) |>
        transmute(contig = X1, length = as.numeric(X2)) |>
        distinct(contig, .keep_all = TRUE) |>
        mutate(sample = s)
})

# ── 3. T2T status: gapless contig with telomeres on both ends ─────────────────
t2t_status <- map_dfr(unique(y_colors$sample), function(s) {
    f <- file.path(QC_DIR, s, "T2T_scaffolds.both.txt")
    if (!file.exists(f)) return(tibble(sample = s, t2t_y = FALSE, scaffold_y = FALSE))
    tbl <- read_tsv(f, col_types = cols(), show_col_types = FALSE) |>
        filter(chromosome == "chrY")
    tibble(
        sample     = s,
        t2t_y      = any(tbl$classification == "contig"),
        scaffold_y = any(tbl$classification == "scaffold")
    )
})

# ── 4. Base layout (all contigs, all samples) ─────────────────────────────────
y_data_all <- y_colors |>
    inner_join(y_lengths, by = c("contig", "sample")) |>
    left_join(t2t_status, by = "sample") |>
    mutate(
        t2t_y      = replace_na(t2t_y, FALSE),
        scaffold_y = replace_na(scaffold_y, FALSE)
    )

# ── Helper: recompute contig positions after subsetting ───────────────────────
build_layout <- function(df, min_length = 0) {
    df |>
        filter(length >= min_length) |>
        group_by(sample) |>
        arrange(desc(length), .by_group = TRUE) |>
        mutate(
            contig_start = cumsum(lag(length, default = 0)),
            total_y      = sum(length),
            n_contigs    = n()
        ) |>
        ungroup()
}

# ── Build the two subsets ─────────────────────────────────────────────────────
y_data_t2t <- y_data_all |>
    filter(t2t_y | scaffold_y) |>
    mutate(bar_class = if_else(t2t_y, "T2T contig", "T2T scaffold")) |>
    build_layout(min_length = 1e6)

y_data_frag <- y_data_all |>
    mutate(bar_class = case_when(
        t2t_y      ~ "T2T contig",
        scaffold_y ~ "T2T scaffold",
        TRUE       ~ "Fragmented"
    )) |>
    build_layout()

# ── Feature loading helpers ───────────────────────────────────────────────────
read_y_bed <- function(s, path, col_names, col_types, active_contigs) {
    if (!file.exists(path)) return(NULL)
    read_tsv(path, col_names = col_names, col_types = col_types,
             show_col_types = FALSE) |>
        filter(contig %in% active_contigs) |>
        mutate(sample = s)
}

join_y_bed <- function(df, layout, sample_ord) {
    layout |>
        select(sample, contig, contig_start) |>
        distinct() |>
        inner_join(df, by = c("sample", "contig")) |>
        mutate(
            abs_start = contig_start + start,
            abs_end   = contig_start + end,
            sample    = factor(sample, levels = sample_ord)
        )
}

load_features <- function(y_data_sub, sample_ord) {
    samps <- unique(as.character(y_data_sub$sample))

    active_contigs <- function(s)
        y_data_sub |> filter(as.character(sample) == s) |> pull(contig)

    dt_gaps <- map_dfr(samps, function(s)
        read_y_bed(s, file.path(QC_DIR, s, "gap_stats.both.n_regions.bed"),
                   c("contig", "start", "end", "label", "length"), "cddcd",
                   active_contigs(s))
    ) |> join_y_bed(y_data_sub, sample_ord)

    dt_satellites <- map_dfr(samps, function(s)
        read_y_bed(s,
                   file.path(RM_DIR, s, "rm_summary", paste0(s, "_filtered_satellites.bed")),
                   c("contig", "start", "end", "family", "name"), "cddcc",
                   active_contigs(s))
    ) |> join_y_bed(y_data_sub, sample_ord)

    dt_telo <- map_dfr(samps, function(s) {
        f <- file.path(QC_DIR, s, "T2T_contigs.both_motif.csv")
        if (!file.exists(f)) return(NULL)
        ac <- active_contigs(s)
        lengths <- y_data_sub |> filter(as.character(sample) == s) |>
            select(contig, length) |> distinct()
        read_csv(f, col_types = "cddd", show_col_types = FALSE) |>
            separate(id, sep = "_", into = c("contig", "location")) |>
            filter(contig %in% ac) |>
            group_by(contig) |>
            mutate(repeat_number = max(forward_repeat_number, reverse_repeat_number)) |>
            ungroup() |>
            filter(repeat_number > 10) |>
            left_join(lengths, by = "contig") |>
            mutate(
                start = if_else(location == "Start", 0L, as.integer(length - repeat_number * 5)),
                end   = if_else(location == "Start", as.integer(repeat_number * 5), length)
            ) |>
            select(contig, start, end) |>
            mutate(sample = s)
    }) |> join_y_bed(y_data_sub, sample_ord)

    list(gaps = dt_gaps, satellites = dt_satellites, telo = dt_telo)
}

# ── Plot function ─────────────────────────────────────────────────────────────
plot_karyogram <- function(y_data_sub, suffix, fill_values, fill_breaks) {
    sample_ord <- y_data_sub |>
        distinct(sample, total_y) |>
        arrange(total_y) |>
        pull(sample) |>
        as.character()

    y_data_sub <- y_data_sub |>
        mutate(
            sample = factor(as.character(sample), levels = sample_ord),
            x_pos  = as.numeric(factor(as.character(sample), levels = sample_ord))
        )

    feats <- load_features(y_data_sub, sample_ord)

    dt_gaps <- feats$gaps |>
        mutate(
            gap_center = (abs_start + abs_end) / 2,
            gap_half   = pmax((abs_end - abs_start) / 2, MIN_GAP_VIS),
            gap_ymin   = gap_center - gap_half,
            gap_ymax   = gap_center + gap_half,
            x_pos      = as.numeric(sample)
        )

    half <- BAR_WIDTH / 2
    n_main <- length(fill_breaks) - 2  # all breaks except Gap and Satellite

    p <- ggplot() +
        geom_rect(
            data = y_data_sub,
            aes(xmin = x_pos - half, xmax = x_pos + half,
                ymin = contig_start,  ymax = contig_start + length,
                fill = bar_class),
            colour = "grey30", linewidth = 0.5
        ) +
        geom_rect(
            data = feats$satellites |> mutate(x_pos = as.numeric(sample)),
            aes(xmin = x_pos - half, xmax = x_pos + half,
                ymin = abs_start,    ymax = abs_end,
                fill = "Satellite"),
            alpha = 0.65, colour = NA
        ) +
        geom_rect(
            data = feats$telo |> mutate(x_pos = as.numeric(sample)),
            aes(xmin = x_pos - half, xmax = x_pos + half,
                ymin = abs_start,    ymax = abs_end),
            fill = "#1b7837", alpha = 0.9, colour = NA
        ) +
        geom_rect(
            data = y_data_sub,
            aes(xmin = x_pos - half, xmax = x_pos + half,
                ymin = contig_start,  ymax = contig_start + length),
            colour = "grey30", fill = NA, linewidth = 0.5
        ) +
        geom_hline(yintercept = REF_Y_LENGTH,
                   linewidth = 0.7, lty = 2, colour = "steelblue") +
        geom_rect(
            data = dt_gaps,
            aes(xmin = x_pos - half * 0.9, xmax = x_pos + half * 0.9,
                ymin = gap_ymin,            ymax = gap_ymax,
                fill = "Gap"),
            colour = "#7f0000", linewidth = 0.3, alpha = 0.9
        ) +
        annotate("text", x = 0.55, y = REF_Y_LENGTH - 1e6,
                 label = sprintf("CHM13 chrY: %.1f Mb", REF_Y_LENGTH / 1e6),
                 hjust = 1, vjust = 0, colour = "steelblue", size = 2.5) +
        coord_flip() +
        scale_x_continuous(
            breaks = seq_along(sample_ord),
            labels = sample_ord,
            expand = expansion(add = 0.6)
        ) +
        scale_y_continuous(
            labels = scales::unit_format(unit = "Mb", scale = 1e-6),
            expand = expansion(mult = c(0, 0.12))
        ) +
        scale_fill_manual(values = fill_values, breaks = fill_breaks,
                          na.value = "grey75") +
        guides(fill = guide_legend(
            nrow = 1, byrow = TRUE,
            override.aes = list(
                colour = c(rep(NA, n_main), "#7f0000", "grey40"),
                alpha  = c(rep(1,  n_main), 0.9,       0.65)
            )
        )) +
        expand_limits(y = c(-5e5, REF_Y_LENGTH * 1.14)) +
        labs(y = "Assembled chrY length", x = NULL, fill = NULL) +
        theme_classic(base_size = 9) +
        theme(
            legend.position    = "bottom",
            legend.key.size    = unit(3, "mm"),
            axis.text.y        = element_text(size = 8),
            axis.line.y        = element_blank(),
            axis.ticks.y       = element_blank(),
            panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
        )

    ggsave(paste0(output_prefix, ".", suffix, ".pdf"), p, width = 7, height = 6)
    ggsave(paste0(output_prefix, ".", suffix, ".png"), p, width = 7, height = 6, dpi = 300)
    cat("Done:", suffix, "\n")
}

# ── Generate plots ────────────────────────────────────────────────────────────
plot_karyogram(
    y_data_sub  = y_data_t2t,
    suffix      = "Y_karyogram",
    fill_values = c("T2T contig"   = "#92b0d1",
                    "T2T scaffold" = "#dab887",
                    "Gap"          = "#d73027",
                    "Satellite"    = "white"),
    fill_breaks = c("T2T contig", "T2T scaffold", "Gap", "Satellite")
)

plot_karyogram(
    y_data_sub  = y_data_frag,
    suffix      = "Y_karyogram_fragmented",
    fill_values = c("T2T contig"   = "#92b0d1",
                    "T2T scaffold" = "#dab887",
                    "Fragmented"   = "#b2abd2",
                    "Gap"          = "#d73027",
                    "Satellite"    = "white"),
    fill_breaks = c("T2T contig", "T2T scaffold", "Fragmented", "Gap", "Satellite")
)
