#!/usr/bin/env Rscript

MIN_LENGTH = 10
MIN_INSERTION_SIZE = 100
MIN_DELETION_SIZE = 100

# Load required libraries
suppressPackageStartupMessages({
    library(SVbyEye)
    library(ggplot2)
    library(dplyr)
    library(cowplot)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript synteny_plot.R <hap1.paf> <hap2.paf> <output.png> <output.pdf> [target_region] [cdna_hap1.paf] [cdna_hap2.paf] [cdna_ref.paf]")
}

hap1_paf <- args[1]
hap2_paf <- args[2]
output_png <- args[3]
output_pdf <- args[4]
cdna_hap1_paf <- if (length(args) >= 6 && args[5] != "") args[6] else NULL
cdna_hap2_paf <- if (length(args) >= 7 && args[6] != "") args[7] else NULL
cdna_ref_paf <- if (length(args) >= 8 && args[7] != "") args[8] else NULL

# Function to sort chromosomes in natural order
sort_chromosomes <- function(chr_names) {
    # Extract chromosome numbers/names
    chr_order <- data.frame(
        chr = chr_names,
        stringsAsFactors = FALSE
    )
    
    # Parse chromosome names - improved regex
    chr_order$number <- gsub("^chr", "", chr_order$chr, ignore.case = TRUE)
    
    # Create sorting key
    chr_order$sort_key <- sapply(chr_order$number, function(x) {
        if (x %in% c("X", "x")) return(23)
        else if (x %in% c("Y", "y")) return(24)
        else if (x %in% c("M", "MT", "m", "mt")) return(25)
        else {
            # Try to convert to numeric, if it fails return a high number
            num_val <- suppressWarnings(as.numeric(x))
            if (is.na(num_val)) return(26)
            else return(num_val)
        }
    })
    
    # Sort and return chromosome names
    chr_order <- chr_order[order(chr_order$sort_key), ]
    return(chr_order$chr)
}

# Function to filter PAF to main chromosomes only (only filters target/reference sequences)
filter_main_chromosomes <- function(paf_table, label = "") {
    message(label, " - Initial alignments: ", nrow(paf_table))
    message(label, " - Unique query contigs: ", length(unique(paf_table$q.name)))
    message(label, " - Unique target sequences: ", length(unique(paf_table$t.name)))
    message(label, " - Target sequences: ", paste(head(unique(paf_table$t.name), 10), collapse = ", "))
    
    # Keep only alignments to main chromosomes (chr1-22, chrX, chrY, chrM)
    # Remove decoy, unlocalized, unplaced, random, alt, HLA sequences
    # This only filters the TARGET (reference) sequences, not the QUERY (assembly contigs)
    main_chr_pattern <- "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$"
    
    paf_filtered <- paf_table %>%
        filter(grepl(main_chr_pattern, t.name, ignore.case = FALSE))
    
    message(label, " - After chromosome filter: ", nrow(paf_filtered), " alignments")
    message(label, " - Unique query contigs remaining: ", length(unique(paf_filtered$q.name)))
    message(label, " - Unique target sequences remaining: ", length(unique(paf_filtered$t.name)))
    message(label, " - Target sequences remaining: ", paste(unique(paf_filtered$t.name), collapse = ", "))
    
    return(paf_filtered)
}

# Function to parse target region (format: chr:start-end)
parse_target_region <- function(region_str) {
    if (is.null(region_str) || region_str == "") {
        return(NULL)
    }
    
    parts <- strsplit(region_str, ":")[[1]]
    if (length(parts) != 2) {
        stop("Invalid region format. Use chr:start-end")
    }
    
    chrom <- parts[1]
    coords <- strsplit(parts[2], "-")[[1]]
    if (length(coords) != 2) {
        stop("Invalid region format. Use chr:start-end")
    }
    
    return(list(
        chr = chrom,
        start = as.numeric(coords[1]),
        end = as.numeric(coords[2])
    ))
}

# Function to subset PAF by target region
subset_by_region <- function(paf_table, region, label = "") {
    if (is.null(region)) {
        return(paf_table)
    }
    
    region_str <- paste0(region$chr, ":", region$start, "-", region$end)
    message(label, " - Subsetting to region: ", region_str)
    message(label, " - Before region filter: ", nrow(paf_table), " alignments")
    
    paf_subset <- subsetPafAlignments(
        paf.table = paf_table,
        target.region = region_str
    )
    
    message(label, " - After region filter: ", nrow(paf_subset), " alignments")
    message(label, " - Unique query contigs after region filter: ", length(unique(paf_subset$q.name)))
    
    return(paf_subset)
}

# Function to write PAF format file
write_paf <- function(paf_table, filename) {
    # PAF format columns
    paf_output <- paf_table %>%
        select(q.name, q.len, q.start, q.end, strand,
               t.name, t.len, t.start, t.end,
               n.match, aln.len, mapq)
    
    # Add CIGAR if available
    if ("cg" %in% colnames(paf_table)) {
        paf_output$cg <- paste0("cg:Z:", paf_table$cg)
        write.table(paf_output, file = filename, quote = FALSE, 
                   sep = "\t", row.names = FALSE, col.names = FALSE)
    } else {
        write.table(paf_output, file = filename, quote = FALSE, 
                   sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

# Function to create summary table
create_summary_table <- function(paf_table) {
    summary_table <- paf_table %>%
        mutate(
            alignment_length = aln.len,
            query_span = q.end - q.start,
            target_span = t.end - t.start,
            cigar = ifelse("cg" %in% colnames(paf_table), cg, NA)
        ) %>%
        select(
            query_name = q.name,
            query_length = q.len,
            query_start = q.start,
            query_end = q.end,
            strand,
            target_name = t.name,
            target_length = t.len,
            target_start = t.start,
            target_end = t.end,
            num_matches = n.match,
            alignment_length,
            query_span,
            target_span,
            mapq,
            cigar
        )
    
    return(summary_table)
}

# Helper function to convert PAF to GenomicRanges
paf_to_granges <- function(paf_table, prefix = "") {
    if (nrow(paf_table) == 0) {
        return(NULL)
    }
    
    # Add prefix to target names if specified
    target_names <- if (prefix != "") {
        paste0(prefix, paf_table$t.name)
    } else {
        paf_table$t.name
    }
    
    # Convert to GenomicRanges
    gr <- GenomicRanges::GRanges(
        seqnames = target_names,
        ranges = IRanges::IRanges(
            start = paf_table$t.start,
            end = paf_table$t.end
        ),
        strand = paf_table$strand,
        ID = paf_table$q.name  # Gene/transcript name
    )
    
    return(gr)
}

# Function to add gene tracks to AVA plot
add_gene_tracks <- function(plt, cdna_hap1, cdna_hap2, cdna_ref) {
    message("\n========== Adding gene tracks ==========")
    
    # Process Hap1 gene alignments
    if (!is.null(cdna_hap1) && file.exists(cdna_hap1)) {
        message("Reading Hap1 cDNA alignments: ", cdna_hap1)
        paf_cdna_hap1 <- readPaf(paf.file = cdna_hap1, include.paf.tags = FALSE)
        
        if (nrow(paf_cdna_hap1) > 0) {
            # Convert to GenomicRanges with Hap1. prefix
            gene_gr <- paf_to_granges(paf_cdna_hap1, prefix = "Hap1.")
            
            if (!is.null(gene_gr)) {
                # Add to plot
                plt <- addAnnotation(
                    ggplot.obj = plt,
                    annot.gr = gene_gr,
                    coordinate.space = "self",
                    y.label.id = "seqnames",
                    shape = "rectangle",
                    annotation.group = "ID",
                    fill.by = "ID",
                    annotation.label = "Genes.Hap1",
                    annotation.level = 0
                )
                message("  Added ", length(gene_gr), " gene alignments for Hap1")
            }
        }
    }
    
    # Process Hap2 gene alignments
    if (!is.null(cdna_hap2) && file.exists(cdna_hap2)) {
        message("Reading Hap2 cDNA alignments: ", cdna_hap2)
        paf_cdna_hap2 <- readPaf(paf.file = cdna_hap2, include.paf.tags = FALSE)
        
        if (nrow(paf_cdna_hap2) > 0) {
            # Convert to GenomicRanges with Hap2. prefix
            gene_gr <- paf_to_granges(paf_cdna_hap2, prefix = "Hap2.")
            
            if (!is.null(gene_gr)) {
                # Add to plot
                plt <- addAnnotation(
                    ggplot.obj = plt,
                    annot.gr = gene_gr,
                    coordinate.space = "self",
                    y.label.id = "seqnames",
                    shape = "rectangle",
                    annotation.group = "ID",
                    fill.by = "ID",
                    annotation.label = "Genes.Hap2",
                    annotation.level = 0
                )
                message("  Added ", length(gene_gr), " gene alignments for Hap2")
            }
        }
    }
    
    # Process Reference gene alignments
    if (!is.null(cdna_ref) && file.exists(cdna_ref)) {
        message("Reading Reference cDNA alignments: ", cdna_ref)
        paf_cdna_ref <- readPaf(paf.file = cdna_ref, include.paf.tags = FALSE)
        
        if (nrow(paf_cdna_ref) > 0) {
            # Convert to GenomicRanges (no prefix for reference)
            gene_gr <- paf_to_granges(paf_cdna_ref, prefix = "")
            
            if (!is.null(gene_gr)) {
                # Add to plot
                plt <- addAnnotation(
                    ggplot.obj = plt,
                    annot.gr = gene_gr,
                    coordinate.space = "self",
                    y.label.id = "seqnames",
                    shape = "rectangle",
                    annotation.group = "ID",
                    fill.by = "ID",
                    annotation.label = "Genes.Ref",
                    annotation.level = 0
                )
                message("  Added ", length(gene_gr), " gene alignments for Reference")
            }
        }
    }
    
    return(plt)
}

# Read PAF files
message("========== Reading PAF files ==========")
message("Reading haplotype 1 PAF file: ", hap1_paf)
paf_hap1 <- readPaf(paf.file = hap1_paf, include.paf.tags = TRUE, restrict.paf.tags = "cg")

message("\nReading haplotype 2 PAF file: ", hap2_paf)
paf_hap2 <- readPaf(paf.file = hap2_paf, include.paf.tags = TRUE, restrict.paf.tags = "cg")

# Filter to main chromosomes
message("\n========== Filtering to main chromosomes ==========")
paf_hap1 <- filter_main_chromosomes(paf_hap1, "Haplotype 1")
message("")
paf_hap2 <- filter_main_chromosomes(paf_hap2, "Haplotype 2")

# Create whole genome plot
message("\n========== Creating whole genome plot ==========")

# Get unique chromosomes from the reference and sort them properly
chromosomes_raw <- unique(c(paf_hap1$t.name, paf_hap2$t.name))
chromosomes <- sort_chromosomes(chromosomes_raw)
message("Chromosomes to plot (sorted): ", paste(chromosomes, collapse = ", "))

# Define consistent color palette for both plots
color_palette <- c("+" = "#D55E00", "-" = "#0072B2")

# Create whole genome plot for haplotype 1
message("\n========== Creating genome plot for Haplotype 1 ==========")
plt_hap1 <- plotGenome(
    paf.table = paf_hap1,
    chromosomes = chromosomes,
    color.by = "direction",
    chromosome.bar.width = grid::unit(3, "mm"),
    min.query.aligned.bp = 100000
)

plt_hap1 <- plt_hap1 +
    scale_fill_manual(
        name = "Strand",
        values = color_palette,
        labels = c("+" = "Forward", "-" = "Reverse"),
        breaks = c("+", "-")
    ) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "Haplotype 1")

# Create whole genome plot for haplotype 2
message("\n========== Creating genome plot for Haplotype 2 ==========")
plt_hap2 <- plotGenome(
    paf.table = paf_hap2,
    chromosomes = chromosomes,
    color.by = "direction",
    chromosome.bar.width = grid::unit(3, "mm"),
    min.query.aligned.bp = 100000
)

plt_hap2 <- plt_hap2 +
    scale_fill_manual(
        name = "Strand",
        values = color_palette,
        labels = c("+" = "Forward", "-" = "Reverse"),
        breaks = c("+", "-")
    ) +
    theme_bw() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(title = "Haplotype 2")

# Extract legend from one of the plots
message("\n========== Extracting shared legend ==========")
plt_legend_temp <- plt_hap1 +
    theme(legend.position = "top",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 11))

shared_legend <- get_legend(plt_legend_temp)

# Combine both plots horizontally with shared legend on top
message("\n========== Combining genome plots ==========")
plots_row <- plot_grid(
    plt_hap1,
    plt_hap2,
    ncol = 2,
    align = "h",
    axis = "tb",
    rel_widths = c(1, 1)
)

# Add shared legend on top
plt <- plot_grid(
    shared_legend,
    plots_row,
    ncol = 1,
    rel_heights = c(0.05, 1)
)

# Save whole genome plots
message("\nSaving whole genome PNG plot to: ", output_png)
ggsave(
    filename = output_png,
    plot = plt,
    width = 24,
    height = 14,
    dpi = 300,
    units = "in"
)

message("Saving whole genome PDF plot to: ", output_pdf)
ggsave(
    filename = output_pdf,
    plot = plt,
    width = 24,
    height = 14,
    units = "in"
)

message("\n========== Done! ==========")
