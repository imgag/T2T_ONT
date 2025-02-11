library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

#' Read and process MinKNOW report data
#' @param csv_path Path to the CSV file containing parsed reports
#' @return Processed tibble with merged flowcell data
process_report_data <- function(csv_path) {
  # Read CSV and convert times to datetime
  data <- read_csv(csv_path) %>%
    mutate(
      run_start_time = as_datetime(run_start_time),
      run_end_time = as_datetime(run_end_time),
      # Convert elapsed time to hours for plotting
      run_time_hours = elapsed_time_seconds / 3600
    )
  
  return(data)
}

#' Create plots for MinKNOW report data
#' @param data Processed report data
#' @param output_dir Directory to save plots
create_plots <- function(data, output_dir) {
  # Common theme elements
  theme_custom <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # Base plot function
  create_metric_plot <- function(data, metric, ylabel, scale = "continuous") {
    p <- ggplot(data, aes(x = run_start_time, y = !!sym(metric), 
                         color = kit_type)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_x_datetime(
        date_breaks = "1 month", 
        date_labels = "%b %Y",
        expand = expansion(mult = c(0.05, 0.05))  # Add small padding
      ) +
      labs(
        x = "Run Start Time",
        y = ylabel,
        color = "Kit Type"
      ) +
      theme_custom
    
    if (scale == "log10") {
      p <- p + scale_y_log10(labels = comma)
    } else {
      p <- p + scale_y_continuous(labels = comma)
    }
    
    return(p)
  }
  
  # Create individual plots
  p1 <- create_metric_plot(data, "estimated_n50", "N50", "log10")
  p2 <- create_metric_plot(data, "bases_called_pass", "Bases Called (Pass)", "log10")
  p3 <- create_metric_plot(data, "reads_called_pass", "Reads Called (Pass)", "log10")
  p4 <- create_metric_plot(data, "run_time_hours", "Run Time (hours)")
  
  # Combine plots
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "MinKNOW Sequencing Run Metrics",
      subtitle = paste("Data from", format(min(data$run_start_time), "%b %Y"), 
                      "to", format(max(data$run_start_time), "%b %Y")),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Save plots
  ggsave(
    file.path(output_dir, "minknow_metrics.pdf"),
    combined_plot,
    width = 15,
    height = 12
  )
  
  ggsave(
    file.path(output_dir, "minknow_metrics.png"),
    combined_plot,
    width = 15,
    height = 12,
    dpi = 300
  )
}

# Main execution
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript plot_minknow_reports.R <csv_file> <output_dir>
    
Description:
  Create plots from MinKNOW report data
  
Arguments:
  csv_file    Path to the CSV file containing parsed report data
  output_dir  Directory to save the output plots\n")
    quit(status = 1)
  }
  
  csv_path <- args[1]
  output_dir <- args[2]
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process data and create plots
  data <- process_report_data(csv_path)
  create_plots(data, output_dir)
  
  cat("Created plots in", output_dir, "\n")
}