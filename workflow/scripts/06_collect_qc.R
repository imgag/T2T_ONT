# Check if running in test mode
test_mode <- FALSE

# Define log file path and setup logging
if (!test_mode) {
  log_file <- file(snakemake@log[[1]], open = "wt")  # Open connection in text write mode
  sink(log_file, type = "output")
  sink(log_file, type = "message")
}

# Load necessary libraries
library(tidyverse)

# Define test dataset values
if (test_mode) {
  input_paf_stat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_stat.haplotype1.txt")
  input_paf_asmstat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmstat.haplotype1.txt")
  input_paf_asmgene <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmgene.haplotype1.txt")
  input_sex <- list("test_data/phased_verkko/published_chr19_duplex15x/sample_sex.txt")
  input_whatshap_stats <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_stats.tsv")
  input_whatshap_compare <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_compare.tsv")
  input_merqury <- list("test_data/phased_verkko/published_chr19_duplex15x/merqury.qv")
  input_findt2t_alignment <- list("test_data/phased_verkko/published_chr19_duplex15x/find_T2T/T2T_contigs.haplotype1_alignment_T2T.txt")
  input_findt2t_motif <- list("test_data/phased_verkko/published_chr19_duplex15x/find_T2T/T2T_contigs.haplotype1_motif_T2T.txt") 
  output_qc_full <- "test_output/qc_full.tsv"
} else {
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
  output_qc_full <- snakemake@output[["tsv"]]
}

# Helper function to extract common metadata from file path
get_file_metadata <- function(file) {
  list(
    haplotype = str_extract(basename(file), "haplotype\\d"),  # Extract haplotype from filename
    asm_method = str_extract(file, "([^/]+)(?=/[^/]+/[^/]+$)"),  # Extract assembly method from path
    asm_name = str_extract(file, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
  )
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

full_table <- bind_rows(
  bind_rows(parse_paf_stat(input_paf_stat)),
  bind_rows(parse_paf_asmstat(input_paf_asmstat)),
  bind_rows(parse_paf_asmgene(input_paf_asmgene)),
  bind_rows(parse_whatshap_stats(input_whatshap_stats)),
  bind_rows(parse_whatshap_compare(input_whatshap_compare)),
  bind_rows(parse_merqury(input_merqury)),
  bind_rows(parse_findt2t_alignment(input_findt2t_alignment)),
  bind_rows(parse_findt2t_motif(input_findt2t_motif)),
  bind_rows(parse_sex(input_sex))  
) 

full_table %>%
  write_tsv(output_qc_full)

# Close the log file connection at the end
if (!test_mode) {
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}


