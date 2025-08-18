library(tidyverse)
library(lubridate)
library(scales)

#' Read and process MinKNOW report data
#' @param csv_path Path to the CSV file containing parsed reports
#' @param bio_sample_path Path to the biological sample TSV file
#' @return Processed tibble with merged flowcell data
process_report_data <- function(csv_path, bio_sample_path) {
  # Read biological sample information
  bio_samples <- read_tsv(bio_sample_path) %>%
    select(name_external, run_flowcell_id)
  
  # Read CSV and convert times to datetime
  data <- read_csv(csv_path) %>%
    filter(kit_type == "SQK-ULK114") %>%
    # Group by flowcell_id and merge metrics
    group_by(flow_cell_id, kit_type) %>%
    summarise(
      run_start_time = min(run_start_time),  # Use earliest start time
      run_end_time = max(run_end_time),      # Use latest end time
      estimated_n50 = weighted.mean(estimated_n50, reads_called_pass),  # Weight N50 by number of reads
      bases_called_pass = sum(bases_called_pass),
      reads_called_pass = sum(reads_called_pass),
      elapsed_time_seconds = sum(elapsed_time_seconds),
      .groups = "drop"
    ) %>%
    # Join with biological sample information
    left_join(bio_samples, by = c("flow_cell_id" = "run_flowcell_id")) %>%
    mutate(
      run_start_time = as_datetime(run_start_time),
      run_end_time = as_datetime(run_end_time),
      run_time_hours = elapsed_time_seconds / 3600,
      # Set default name for any unmatched flowcells
      name_external = if_else(is.na(name_external), "Unknown", name_external)
    ) %>%
    # Create long format data for faceting
    pivot_longer(
      cols = c(estimated_n50, bases_called_pass),
      names_to = "metric",
      values_to = "value"
    ) %>%
    # Add nice labels for facets
    mutate(
      metric = factor(metric,
        levels = c("estimated_n50", "bases_called_pass"),
        labels = c("N50", "Bases Called (Pass)")
      )
    )
  
  return(data)
}

#' Create plots for MinKNOW report data
#' @param data Processed report data
#' @param output_dir Directory to save plots
create_plots <- function(data, output_dir) {
  # Create faceted plot
  p <- ggplot(data, aes(x = run_start_time, y = value, 
                        fill = name_external, 
                        shape = kit_type)) +
    geom_point(alpha = 0.7, size = 3) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    # Add both month and week ticks
    scale_x_datetime(
      # Major breaks for months with labels
      date_breaks = "1 month",
      date_labels = "%b %Y",
      # Minor breaks for weeks (ticks only)
      minor_breaks = "1 week",
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(labels = comma) +
    labs(
      x = "Run Start Time",
      y = NULL,
      fill = "Sample",
      shape = "Kit Type",
      title = "MinKNOW Sequencing Run Metrics",
      subtitle = paste("Data from", format(min(data$run_start_time), "%b %Y"),
                      "to", format(max(data$run_start_time), "%b %Y"))
    ) +
    theme_classic() +
    scale_fill_manual("Sample", values = rep(c(
      "#E69F00", "#56B4E9", "#009E73", 
      "#F0E442", "#0072B2", "#D55E00", 
      "#CC79A7", "#999999",
      "#882255", "#44AA99", "#117733",
      "#332288", "#88CCEE", "#DDCC77",
      "#AA4499", "#6699CC"
    ), 2)) +
    scale_shape_manual(values = c(
      "SQK-ULK114" = 21,  # Filled circle
      "SQK-LSK114" = 22,  # Filled triangle
      "SQK-LSK114-XL" = 23,  # Filled square 
      "SQK-APK114" = 24   # Filled diamond
    )) +
    guides(fill = guide_legend("Sample", override.aes = list(shape = 21))) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      # Remove all gridlines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Show axis ticks for weeks
      axis.ticks.x.minor = element_line(size = 0.25),
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Save plots
  ggsave(
    file.path(output_dir, "minknow_metrics_UL.pdf"),
    p,
    width = 12,
    height = 8
  )
  
  ggsave(
    file.path(output_dir, "minknow_metrics_UL.png"),
    p,
    width = 12,
    height = 6,
    dpi = 300
  )
}

# Main execution
if (sys.nframe() == 0) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 3) {
    cat("Usage: Rscript plot_minknow_reports.R <csv_file> <bio_sample_file> <output_dir>
    
Description:
  Create plots from MinKNOW report data
  
Arguments:
  csv_file         Path to the CSV file containing parsed report data
  bio_sample_file  Path to the TSV file containing biological sample information
  output_dir       Directory to save the output plots\n")
    quit(status = 1)
  }
  
  csv_path <- args[1]
  bio_sample_path <- args[2]
  output_dir <- args[3]
  
  # Add test_Data
  csv_path <- "doc/run_reports/run_summary.csv"
  bio_sample_path <- "doc/tables/flowcell_biological_sample.tsv"
  output_dir <- "doc/img"

  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process data and create plots
  data <- process_report_data(csv_path, bio_sample_path)
  create_plots(data, output_dir)
  
  cat("Created plots in", output_dir, "\n")
}