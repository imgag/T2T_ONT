#!/usr/bin/env Rscript
"""
Parse lifted GFF3 files from Flagger/Nucflag annotations and create a merged TSV for plotting.
"""

library(optparse)
library(dplyr)
library(readr)
library(stringr)

# Command line options
option_list <- list(
  make_option(c("--input_dir"), type="character", default="analysis_other/annotations",
              help="Directory containing subdirectories with lifted GFF3 files [default: %default]"),
  make_option(c("--pattern"), type="character", default="flagger_nucflag_annotations.lifted.gff3",
              help="Pattern to match GFF3 files [default: %default]"),
  make_option(c("--output"), type="character", default="analysis_other/annotations/merged_difficult_regions.tsv",
              help="Output TSV file [default: %default]"),
  make_option(c("--samples"), type="character", default=NULL,
              help="Comma-separated list of samples to include (default: all found)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to parse GFF3 attributes
parse_gff3_attributes <- function(attributes_string) {
  # Split by semicolon and parse key=value pairs
  pairs <- str_split(attributes_string, ";")[[1]]
  result <- list()
  
  for (pair in pairs) {
    if (str_detect(pair, "=")) {
      kv <- str_split(pair, "=", n=2)[[1]]
      if (length(kv) == 2) {
        key <- str_trim(kv[1])
        value <- str_trim(kv[2])
        result[[key]] <- value
      }
    }
  }
  return(result)
}

# Function to read and parse a single GFF3 file
read_gff3_file <- function(filepath, sample_name) {
  cat("Processing:", filepath, "\n")
  
  # Read the file, skipping comment lines
  lines <- readLines(filepath)
  data_lines <- lines[!str_starts(lines, "##") & !str_starts(lines, "#")]
  
  if (length(data_lines) == 0) {
    cat("Warning: No data lines found in", filepath, "\n")
    return(data.frame())
  }
  
  # Parse GFF3 format: seqname, source, feature, start, end, score, strand, frame, attributes
  gff_data <- read_tsv(paste(data_lines, collapse="\n"), 
                       col_names = c("seqname", "source", "feature", "start", "end", 
                                   "score", "strand", "frame", "attributes"),
                       col_types = "ccciicccc",
                       show_col_types = FALSE)
  
  if (nrow(gff_data) == 0) {
    return(data.frame())
  }
  
  # Parse attributes for each row
  parsed_data <- gff_data %>%
    rowwise() %>%
    mutate(
      sample = sample_name,
      # Parse key attributes
      attrs = list(parse_gff3_attributes(attributes)),
      # Extract common attributes
      ID = ifelse("ID" %in% names(attrs), attrs[["ID"]], NA_character_),
      Name = ifelse("Name" %in% names(attrs), attrs[["Name"]], NA_character_),
      feature_type = ifelse("feature_type" %in% names(attrs), attrs[["feature_type"]], feature),
      description = ifelse("description" %in% names(attrs), attrs[["description"]], NA_character_),
      coverage = ifelse("coverage" %in% names(attrs), as.numeric(attrs[["coverage"]]), NA_real_),
      sequence_ID = ifelse("sequence_ID" %in% names(attrs), as.numeric(attrs[["sequence_ID"]]), NA_real_),
      extra_copy_number = ifelse("extra_copy_number" %in% names(attrs), as.numeric(attrs[["extra_copy_number"]]), NA_real_),
      low_identity = ifelse("low_identity" %in% names(attrs), attrs[["low_identity"]] == "True", FALSE),
      partial_mapping = ifelse("partial_mapping" %in% names(attrs), attrs[["partial_mapping"]] == "True", FALSE),
      color = ifelse("color" %in% names(attrs), attrs[["color"]], NA_character_),
      # Calculate region size
      region_size = end - start + 1,
      # Standardize feature types
      feature_category = case_when(
        str_detect(tolower(feature_type), "err") ~ "Error",
        str_detect(tolower(feature_type), "dup") ~ "Duplication", 
        str_detect(tolower(feature_type), "hap") ~ "Haploid",
        str_detect(tolower(feature_type), "col") ~ "Collapsed",
        str_detect(tolower(feature_type), "good") ~ "Good",
        str_detect(tolower(feature), "gap") ~ "Gap",
        TRUE ~ "Other"
      ),
      # Assign colors based on feature category if not provided
      plot_color = case_when(
        !is.na(color) ~ paste0("#", str_remove(color, "^#?")),
        feature_category == "Error" ~ "#A20025",      # Dark red
        feature_category == "Duplication" ~ "#FA6800", # Orange  
        feature_category == "Haploid" ~ "#008A00",     # Green
        feature_category == "Collapsed" ~ "#AA00FF",   # Purple
        feature_category == "Good" ~ "#00AA00",        # Light green
        feature_category == "Gap" ~ "#808080",         # Gray
        TRUE ~ "#CCCCCC"                               # Light gray
      )
    ) %>%
    ungroup() %>%
    select(-attrs, -attributes)
  
  return(as.data.frame(parsed_data))
}

# Main execution
main <- function() {
  cat("Starting GFF3 parsing...\n")
  cat("Input directory:", opt$input_dir, "\n")
  cat("Pattern:", opt$pattern, "\n")
  
  # Find all GFF3 files
  if (!is.null(opt$samples)) {
    samples <- str_trim(str_split(opt$samples, ",")[[1]])
    gff_files <- file.path(opt$input_dir, samples, opt$pattern)
    gff_files <- gff_files[file.exists(gff_files)]
    sample_names <- samples[file.exists(file.path(opt$input_dir, samples, opt$pattern))]
  } else {
    # Auto-discover samples
    subdirs <- list.dirs(opt$input_dir, full.names = FALSE, recursive = FALSE)
    gff_files <- file.path(opt$input_dir, subdirs, opt$pattern)
    existing_files <- file.exists(gff_files)
    gff_files <- gff_files[existing_files]
    sample_names <- subdirs[existing_files]
  }
  
  if (length(gff_files) == 0) {
    stop("No GFF3 files found matching pattern: ", opt$pattern)
  }
  
  cat("Found", length(gff_files), "GFF3 files:\n")
  for (i in seq_along(gff_files)) {
    cat("  ", sample_names[i], ":", gff_files[i], "\n")
  }
  
  # Process each file
  all_data <- list()
  for (i in seq_along(gff_files)) {
    tryCatch({
      data <- read_gff3_file(gff_files[i], sample_names[i])
      if (nrow(data) > 0) {
        all_data[[sample_names[i]]] <- data
      }
    }, error = function(e) {
      cat("Error processing", gff_files[i], ":", e$message, "\n")
    })
  }
  
  if (length(all_data) == 0) {
    stop("No data could be parsed from any GFF3 files")
  }
  
  # Combine all data
  merged_data <- bind_rows(all_data)
  
  # Create summary statistics
  cat("\n=== SUMMARY ===\n")
  cat("Total regions:", nrow(merged_data), "\n")
  cat("Samples:", length(unique(merged_data$sample)), "\n")
  cat("Feature categories:\n")
  print(table(merged_data$feature_category))
  
  cat("\nChromosome distribution:\n")
  chrom_summary <- merged_data %>%
    count(seqname, sort = TRUE) %>%
    head(10)
  print(chrom_summary)
  
  cat("\nRegion size statistics (bp):\n")
  size_stats <- merged_data %>%
    group_by(feature_category) %>%
    summarise(
      count = n(),
      mean_size = round(mean(region_size)),
      median_size = round(median(region_size)),
      total_bp = sum(region_size),
      .groups = "drop"
    )
  print(size_stats)
  
  # Prepare final output columns
  output_data <- merged_data %>%
    select(
      sample, seqname, source, feature, start, end, strand,
      feature_type, feature_category, description,
      region_size, coverage, sequence_ID, extra_copy_number,
      low_identity, partial_mapping, plot_color,
      ID, Name
    ) %>%
    arrange(sample, seqname, start)
  
  # Create output directory if needed
  output_dir <- dirname(opt$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write output
  write_tsv(output_data, opt$output)
  cat("\nMerged data written to:", opt$output, "\n")
  cat("Output contains", nrow(output_data), "regions from", 
      length(unique(output_data$sample)), "samples\n")
  
  # Write a summary file as well
  summary_file <- str_replace(opt$output, "\\.tsv$", "_summary.tsv")
  write_tsv(size_stats, summary_file)
  cat("Summary statistics written to:", summary_file, "\n")
}

# Run main function
if (!interactive()) {
  main()
}