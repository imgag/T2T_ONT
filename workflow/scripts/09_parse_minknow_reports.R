library(tidyverse)
library(jsonlite)

#' Extract JSON data from HTML file
#' @param file_path Path to HTML file
#' @return JSON string containing report data
extract_json_from_html <- function(file_path) {
  # Read HTML file
  html_content <- readLines(file_path)
  
  # Find line with reportData
  json_line <- html_content[grep("const reportData=", html_content)]
  
  # Extract JSON string
  json_str <- str_extract(json_line, "(?<=const reportData=).*")
  
  return(json_str)
}

#' Parse report data into tibbles
#' @param file_path Path to HTML file
#' @return List of tibbles containing parsed data

#' Parse report data into tibbles
#' @param file_path Path to HTML file
#' @return List of tibbles containing parsed data
parse_report <- function(file_path) {
  # Extract JSON
  json_str <- extract_json_from_html(file_path)
  
  if (is.null(json_str)) {
    warning(sprintf("Skipping %s due to extraction error", file_path))
    return(NULL)
  }
  
  # Parse JSON with error handling
  report_data <- tryCatch({
    fromJSON(json_str)
  }, error = function(e) {
    warning(sprintf("JSON parsing error in %s: %s", file_path, e$message))
    return(NULL)
  })
  
  if (is.null(report_data)) {
    return(NULL)
  }
  
  # Create summary tibble with all single-value metrics
  run_summary <- tibble(
    # Run info
    flow_cell_id = report_data$flow_cell_id,
    run_start_time = as.POSIXct(report_data$run_start_time),
    run_end_time = as.POSIXct(report_data$run_end_time),
    run_duration_hours = as.numeric(difftime(run_end_time, run_start_time, units = "hours")),
    run_status = report_data$run_status,
    run_status_context = report_data$run_status_additional_context,
    run_complete = report_data$run_complete,
    estimated_n50 = report_data$estimated_n50,
    
    # Device info
    device_type = report_data$header$device_type,
    device_serial = report_data$header$serial,
    experiment_name = report_data$header$experiment_name,
    sample_id = report_data$header$sample_id,
    position = report_data$header$position,
    protocol_run_id = report_data$header$protocol_run_id,
    
    # Feature flags
    alignment_enabled = report_data$alignment_enabled,
    barcoding_enabled = report_data$barcoding_enabled,
    basecalling_enabled = report_data$basecalling_enabled,
    duplex_enabled = report_data$duplex_enabled,
    
    # Output stats
    estimated_bases = report_data$data_output$estimated_bases,
    data_produced_bytes = report_data$data_output$data_produced,
    reads_generated = report_data$data_output$reads_generated,
    
    # Basecalling stats
    reads_called_percent = report_data$basecalling$reads_called,
    reads_called_pass = report_data$basecalling$reads_called_pass,
    reads_called_fail = report_data$basecalling$reads_called_fail,
    reads_called_skipped = report_data$basecalling$reads_called_skipped,
    bases_called_pass = report_data$basecalling$bases_called_pass,
    
    # Run limits
    target_bases = report_data$run_until$target_estimated_bases,
    target_basecalled = report_data$run_until$target_basecalled_bases,
    pores_remaining = report_data$run_until$pores_remaining,
    run_time_limit_seconds = report_data$run_until$total_run_time,
    elapsed_time_seconds = report_data$run_until$elapsed_time_since_start_seconds
  )
  
  # Create time series tibble for metrics over time
  time_series <- tibble(
    flow_cell_id = report_data$flow_cell_id,
    experiment_name = report_data$header$experiment_name,
    sample_id = report_data$header$sample_id,
    
    # Add any time series data here if available in the JSON
    # Example:
    # time = report_data$time_series$time,
    # value = report_data$time_series$value,
    # metric = "example_metric"
  )
  
  return(list(
    summary = run_summary,
    time_series = time_series
  ))
}

#' Parse multiple report files
#' @param report_dir Directory containing report HTML files
#' @return List of combined tibbles
parse_reports <- function(report_dir) {
  # Find all HTML files
  report_files <- list.files(report_dir, pattern = "*.html", full.names = TRUE)
  
  if (length(report_files) == 0) {
    stop("No HTML files found in ", report_dir)
  }
  
  # Parse each file with error handling
  reports <- map(report_files, parse_report)
  
  # Remove NULL results from failed parses
  reports <- compact(reports)
  
  if (length(reports) == 0) {
    stop("No valid reports could be parsed")
  }
  
  # Combine tibbles
  combined <- list(
    summary = bind_rows(map(reports, "summary")),
    time_series = bind_rows(map(reports, "time_series"))
  )
  
  return(combined)
}

# Example usage:
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("Usage: Rscript parse_reports.R <report_directory>
    
Description:
  Parse MinKNOW report HTML files and extract data into CSV files.
  
Arguments:
  report_directory    Directory containing MinKNOW HTML report files

Output:
  Creates two CSV files in the report directory:
  - run_summary.csv:  One row per sequencing run with all metrics
  - time_series.csv:  Long format time series data\n")
    quit(status = 1)
  }
  
  report_dir <- args[1]
  
  reports <- parse_reports(report_dir)
  
  # Save tibbles as CSV
  write_csv(reports$summary, file.path(report_dir, "run_summary.csv"))
  write_csv(reports$time_series, file.path(report_dir, "time_series.csv"))
}