# Check if running in test mode
test_mode <- FALSE

# Define log file path and setup logging
if (!test_mode) {
  log_file <- file(snakemake@log[[1]], open = "wt")  # Open connection in text write mode
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

# Start time logging
start_time <- Sys.time()
message(paste("Starting QC collection at:", start_time))

# Load necessary libraries
library(tidyverse)
message("Libraries loaded successfully")

# Define test dataset values
if (test_mode) {
  message("Running in test mode")
  input_paf_stat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_stat.haplotype1.txt")
  input_paf_asmstat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmstat.haplotype1.txt")
  input_paf_asmgene <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmgene.haplotype1.txt")
  input_sex <- list("test_data/phased_verkko/published_chr19_duplex15x/sample_sex.txt")
  input_whatshap_stats <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_stats.tsv")
  input_whatshap_compare <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_compare.tsv")
  input_merqury <- list("test_data/phased_verkko/published_chr19_duplex15x/merqury.qv")
  input_findt2t_alignment <- list("test_data/phased_verkko/published_chr19_duplex15x/find_T2T/T2T_contigs.haplotype1_alignment_T2T.txt")
  input_findt2t_motif <- list("test_data/phased_verkko/published_chr19_duplex15x/find_T2T/T2T_contigs.haplotype1_motif_T2T.txt") 
  input_gap_stats <- list("test_data/phased_verkko/published_chr19_duplex15x/n_stats.haplotype1.tsv")
  output_qc_full <- "test_output/qc_full.tsv"
} else {
  message("Running in production mode")
  # Access Snakemake input and output
  input_paf_stat <- snakemake@input[["paf_stat"]]
  input_paf_asmstat <- snakemake@input[["paf_asmstat"]]
  input_paf_asmgene <- snakemake@input[["paf_asmgene"]]
  input_sex <- snakemake@input[["sex"]]
  input_whatshap_stats <- snakemake@input[["whatshap_stats"]]
  input_whatshap_compare <- snakemake@input[["whatshap_compare"]]
  input_merqury <- snakemake@input[["merqury_stats"]]
  input_findt2t_alignment <- snakemake@input[["findt2t_alignment"]]
  input_findt2t_motif <- snakemake@input[["findt2t_motif"]]
  input_gap_stats <- snakemake@input[["gap_stats"]]
  input_contig_stats <- snakemake@input[["contig_stats"]]
  output_qc_full <- snakemake@output[["tsv"]]
}

message("Input files loaded:")
message(paste("- PAF stats:", length(input_paf_stat), "files"))
message(paste("- PAF asmstat:", length(input_paf_asmstat), "files"))
message(paste("- PAF asmgene:", length(input_paf_asmgene), "files"))
message(paste("- Sex determination:", length(input_sex), "files"))
message(paste("- WhatsHap stats:", length(input_whatshap_stats), "files"))
message(paste("- WhatsHap compare:", length(input_whatshap_compare), "files"))
message(paste("- Merqury:", length(input_merqury), "files"))
message(paste("- FindT2T alignment:", length(input_findt2t_alignment), "files"))
message(paste("- FindT2T motif:", length(input_findt2t_motif), "files"))
message(paste("- Gap statistics:", length(input_gap_stats), "files"))
message(paste("- Contig statistics:", length(input_contig_stats), "files"))
# Helper function to extract common metadata from file path
get_file_metadata <- function(file) {
  list(
    haplotype = str_extract(basename(file), "haplotype\\d"),  # Extract haplotype from filename
    asm_method = str_extract(file, "([^/]+)(?=/[^/]+/[^/]+$)"),  # Extract assembly method from path
    asm_name = str_extract(file, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
  )
}

# Function to time and log execution of parsing functions
time_execution <- function(func, args, func_name) {
  message(paste("Starting", func_name, "at:", Sys.time()))
  start <- Sys.time()
  result <- func(args)
  end <- Sys.time()
  time_taken <- end - start
  message(paste("Completed", func_name, "in:", round(time_taken, 2), "seconds"))
  return(result)
}

# Example of simplified parsing function
parse_paf_stat <- function(paf_stats) {
  results <- lapply(paf_stats, function(file) {
    metadata <- get_file_metadata(file)
    
    read_lines(file) %>%
      enframe(name = NULL, value = "line") %>%
      separate(line, into = c("metric", "value"), sep = ": ", convert = TRUE) %>%
      mutate(metric = str_trim(metric),
             value = as.numeric(value),
             haplotype = metadata$haplotype,
             asm_method = metadata$asm_method,
             asm_name = metadata$asm_name,
             source = "stat")
  })
  names(results) <- paf_stats
  return(results)
}

parse_paf_asmstat <- function(paf_stats) {
  print(paste("Parse paf_asmstat from ", paf_stats))
  # Use lapply to read each file and convert it into a tibble
  results <- lapply(paf_stats, function(paf_stat) {
    metadata <- get_file_metadata(paf_stat)
    
    # Read the file as a table
    read_tsv(paf_stat, col_names = c("metric", "value"), skip = 1, show_col_types = FALSE) %>%
      mutate(metric = str_trim(metric), 
             value = case_when(         
               grepl("%$", value) ~ as.numeric(sub("%", "", value)) / 100,
               TRUE ~ suppressWarnings(as.numeric(value))
             ),
             haplotype = metadata$haplotype,
             asm_method = metadata$asm_method,
             asm_name = metadata$asm_name,
             source = "asmstat")
  })
  
  # Set names for the results list based on the input filenames
  names(results) <- paf_stats
  return(results)
}

parse_paf_asmgene <- function(paf_stats) {
  print(paste("Parse paf_asmgene from ", paf_stats))
  # Use lapply to read each file and convert it into a tibble
  results <- lapply(paf_stats, function(paf_stat) {
    metadata <- get_file_metadata(paf_stat)
    
    # Read the file as a table and transform to long format
    read_tsv(paf_stat, col_names = c("col_type", "metric", "ref", "asm"), skip = 1, show_col_types = FALSE) %>%
      select(-col_type) %>%
      pivot_longer(cols = c(ref, asm), names_to = "col", values_to = "value") %>%
      mutate(metric = str_c(col, ".", str_trim(metric)),  # Add prefix to metric
             value = case_when(         
               grepl("%$", value) ~ as.numeric(sub("%", "", value)) / 100,
               TRUE ~ suppressWarnings(as.numeric(value))
             ),
             haplotype = metadata$haplotype,
             asm_method = metadata$asm_method,
             asm_name = metadata$asm_name,
             source = "asmgene") %>%
      select(-col)
  })
  
  # Set names for the results list based on the input filenames
  names(results) <- paf_stats
  return(results)
}

parse_sex <- function(files) {
  print(paste("Parse sex from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    
    # Read single value from file
    sex_value <- read_lines(file, n_max = 1)
    sex_value <- str_split_1(sex_value, "\t")[1]
    # Convert to numeric (0 = female, 1 = male)
    numeric_value <- case_when(
      sex_value == "female" ~ 0,
      sex_value == "male" ~ 1,
      TRUE ~ NA_real_  # For "unknown" or any other value
    )
    
    tibble(
      metric = "n_y_chrom",
      value = numeric_value,
      haplotype = "both",
      asm_method = metadata$asm_method,
      asm_name = metadata$asm_name,
      source = "sex_determination"
    )
  })
  names(results) <- files
  return(results)
}

parse_whatshap_stats <- function(file) {
  print(paste("Parse whatshap_stats from ", file))
  results <- lapply(file, function(file) {
    metadata <- get_file_metadata(file)
    read_tsv(file, col_names = TRUE, show_col_types = FALSE) %>%
      select(c(-`#sample`, -file_name)) %>%
      pivot_longer(cols = c(-chromosome), names_to = "metric", values_to = "value") %>%
      mutate(metric = metric,
             value = as.numeric(value),
             haplotype = "both",
             asm_method = metadata$asm_method,
             asm_name = metadata$asm_name,
             source = "whatshap_stats")
  })
  names(results) <- file
  return(results)
}

parse_whatshap_compare <- function(file) {
  print(paste("Parse whatshap_compare from ", file))
  results <- lapply(file, function(file) {
    metadata <- get_file_metadata(file)
    read_tsv(file, col_names = TRUE, show_col_types = FALSE) %>%
      select(c(-`#sample`, -dataset_name0, -dataset_name1, -file_name0, -file_name1, -all_switchflips, -largestblock_switchflips)) %>%
      pivot_longer(cols = c(-chromosome), names_to = "metric", values_to = "value") %>%
      mutate(metric = metric,
             value = as.numeric(value),
             haplotype = "both",
             asm_method = metadata$asm_method,
             asm_name = metadata$asm_name,
             source = "whatshap_compare")
  })
  names(results) <- file
  return(results)
}

parse_merqury <- function(files) {
  print(paste("Parse merqury from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    
    read_tsv(file, col_names = c("haplotype_name", "unique_kmers", "total_kmers", "qv", "error_rate"), show_col_types = FALSE) %>%
      pivot_longer(cols = c("unique_kmers", "total_kmers", "qv", "error_rate"), 
                  names_to = "metric", 
                  values_to = "value") %>%
      mutate(
        metric = str_c(metric),
        haplotype = haplotype_name,  # Use the actual haplotype from the file
        asm_method = metadata$asm_method,
        asm_name = metadata$asm_name,
        source = "merqury"
      ) %>%
      select(-haplotype_name)  # Remove the temporary column
  })
  names(results) <- files
  return(results)
}

parse_findt2t_alignment <- function(files) {
  print(paste("Parse findt2t_alignment from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    
    # Count number of lines (alignments)
    n_t2t <- nrow(read_tsv(
      file, show_col_types = FALSE,
      col_names = FALSE))
    
    # New metric: Count number of T2T chromosomes
    tibble(
      metric = "n_T2T",
      value = n_t2t,
      haplotype = metadata$haplotype,
      asm_method = metadata$asm_method,
      asm_name = metadata$asm_name,
      source = "t2t_alignment"
    )
  })
  names(results) <- files
  return(results)
}

parse_findt2t_motif <- function(files) {
  print(paste("Parse findt2t_motif from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    
    # Count number of lines (motifs)
    n_motif <- nrow(read_tsv(file, show_col_types = FALSE))
    
    tibble(
      metric = "n_T2T",
      value = n_motif,
      haplotype = metadata$haplotype,
      asm_method = metadata$asm_method,
      asm_name = metadata$asm_name,
      source = "t2t_motif"
    )
  })
  names(results) <- files
  return(results)
}

parse_gap_stats <- function(files) {
  print(paste("Parse N statistics from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    
    read_tsv(file, col_names = c("metric", "value"), show_col_types = FALSE) %>%
      # Skip header row if it exists
      filter(metric != "Metric") %>%
      mutate(
        # Clean up metric names for consistency
        metric = case_when(
          metric == "Total N count" ~ "total_n_count",
          metric == "Total gaps (N regions)" ~ "total_gaps",
          metric == "Average gap size" ~ "avg_gap_size",
          metric == "Gap N50" ~ "gap_n50",
          metric == "Maximum gap size" ~ "max_gap_size",
          metric == "Minimum gap size" ~ "min_gap_size",
          TRUE ~ metric
        ),
        value = as.numeric(value),
        haplotype = metadata$haplotype,
        asm_method = metadata$asm_method,
        asm_name = metadata$asm_name,
        source = "gap_statistics"
      )
  })
  names(results) <- files
  return(results)
}

parse_contig_stats <- function(files) {
  print(paste("Parse N statistics from ", files))
  results <- lapply(files, function(file) {
    metadata <- get_file_metadata(file)
    read_tsv(file, col_names = c("metric", "value"), show_col_types = FALSE) %>%
      mutate(
        value = as.numeric(value),
        haplotype = metadata$haplotype,
        asm_method = metadata$asm_method,
        asm_name = metadata$asm_name,
        source = "contig_statistics"
      )
  })
  names(results) <- files
  return(results)
}

message("Starting data parsing at:", Sys.time())

# Use the time_execution function for each parsing step
parsed_paf_stat <- time_execution(parse_paf_stat, input_paf_stat, "PAF stat parsing")
parsed_paf_asmstat <- time_execution(parse_paf_asmstat, input_paf_asmstat, "PAF asmstat parsing")
parsed_paf_asmgene <- time_execution(parse_paf_asmgene, input_paf_asmgene, "PAF asmgene parsing")
parsed_whatshap_stats <- time_execution(parse_whatshap_stats, input_whatshap_stats, "WhatsHap stats parsing")
parsed_whatshap_compare <- time_execution(parse_whatshap_compare, input_whatshap_compare, "WhatsHap compare parsing")
parsed_merqury <- time_execution(parse_merqury, input_merqury, "Merqury parsing")
parsed_findt2t_alignment <- time_execution(parse_findt2t_alignment, input_findt2t_alignment, "FindT2T alignment parsing")
parsed_findt2t_motif <- time_execution(parse_findt2t_motif, input_findt2t_motif, "FindT2T motif parsing")
parsed_sex <- time_execution(parse_sex, input_sex, "Sex determination parsing")
parsed_gap_stats <- time_execution(parse_gap_stats, input_gap_stats, "Gap statistics parsing")
parsed_contig_stats <- time_execution(parse_contig_stats, input_contig_stats, "Contig statistics parsing")

message("Starting data binding at:", Sys.time())
binding_start <- Sys.time()

full_table <- bind_rows(
  bind_rows(parsed_paf_stat),
  bind_rows(parsed_paf_asmstat),
  bind_rows(parsed_paf_asmgene),
  bind_rows(parsed_whatshap_stats),
  bind_rows(parsed_whatshap_compare),
  bind_rows(parsed_merqury),
  bind_rows(parsed_findt2t_alignment),
  bind_rows(parsed_findt2t_motif),
  bind_rows(parsed_sex),
  bind_rows(parsed_gap_stats),
  bind_rows(parsed_contig_stats) 
) 

binding_end <- Sys.time()
message(paste("Data binding completed in:", round(binding_end - binding_start, 2), "seconds"))

# Write_output
message("Writing output file at:", Sys.time())
write_start <- Sys.time()

full_table %>%
  write_tsv(output_qc_full, col_names = TRUE)

write_end <- Sys.time()
message(paste("Output file written in:", round(write_end - write_start, 2), "seconds"))

# Calculate and report total execution time
end_time <- Sys.time()
total_time <- end_time - start_time
message(paste("Total execution time:", round(total_time, 2), "seconds"))
message(paste("QC collection completed at:", end_time))

# Summary statistics
message("Summary statistics:")
message(paste("- Total metrics collected:", nrow(full_table)))
message(paste("- Unique metric types:", length(unique(full_table$metric))))
message(paste("- Number of haplotypes:", length(unique(full_table$haplotype))))
message(paste("- Number of data sources:", length(unique(full_table$source))))

# Close the log file connection at the end
if (!test_mode) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}


