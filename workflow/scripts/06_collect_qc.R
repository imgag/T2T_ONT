# Load necessary libraries
library(tidyverse)

# Check if running in test mode
test_mode <- TRUE

# Define log file path
if (!test_mode) {
  log_file <- snakemake@log[[1]]
}

# Redirect stdout and stderr to the log file
if (!test_mode) {
  sink(log_file, append = TRUE, type = "output")
  sink(log_file, append = TRUE, type = "message")
}

# Define test dataset values
if (test_mode) {
  input_paf_stat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_stat.haplotype1.txt")
  input_paf_asmstat <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmstat.haplotype1.txt")
  input_paf_asmgene <- list("test_data/phased_verkko/published_chr19_duplex15x/qc_paftools_asmgene.haplotype1.txt")
  input_sex <- list("test_data/phased_verkko/published_chr19_duplex15x/sample_sex.txt")
  input_whatshap_stats <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_stats.tsv")
  input_whatshap_compare <- list("test_data/phased_verkko/published_chr19_duplex15x/whatshap_compare.tsv")
  output_qc_full <- "test_output/qc_full.tsv"
} else {
  # Access Snakemake input and output
  input_paf_stat <- snakemake@input[["paf_stat"]]
  input_paf_asmstat <- snakemake@input[["paf_asmstat"]]
  input_paf_asmgene <- snakemake@input[["paf_asmgene"]]
  input_sex <- snakemake@input[["sex"]]
  input_whatshap_stats <- snakemake@input[["whatshap_stats"]]
  input_whatshap_compare <- snakemake@input[["whatshap_compare"]]
  output_qc_full <- snakemake@output[["qc_full"]]
}


parse_paf_stat <- function(paf_stats) {
  # Use lapply to read each file and convert it into a tibble
  results <- lapply(paf_stats, function(paf_stat) {
    haplotype <- str_extract(basename(paf_stat), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(paf_stat, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(paf_stat, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    read_lines(paf_stat) %>%
      enframe(name = NULL, value = "line") %>%
      separate(line, into = c("metric", "value"), sep = ": ", convert = TRUE) %>%
      mutate(metric = str_trim(metric),
             value = as.numeric(value),
             haplotype = haplotype,
             asm_method = asm_method,
             asm_name = asm_name,
             source = "stat")
  })
  
  # Set names for the results list based on the input filenames
  names(results) <- paf_stats
  return(results)
}


parse_paf_asmstat <- function(paf_stats) {
  # Use lapply to read each file and convert it into a tibble
  results <- lapply(paf_stats, function(paf_stat) {
    haplotype <- str_extract(basename(paf_stat), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(paf_stat, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(paf_stat, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    
    # Read the file as a table
    read_tsv(paf_stat, col_names = c("metric", "value"), skip = 1) %>%
      mutate(metric = str_trim(metric), 
             value = case_when(         
               grepl("%$", value) ~ as.numeric(sub("%", "", value)) / 100,
               TRUE ~ suppressWarnings(as.numeric(value))
             ),
             haplotype = haplotype,
             asm_method = asm_method,
             asm_name = asm_name,
             source = "asmstat")
  })
  
  # Set names for the results list based on the input filenames
  names(results) <- paf_stats
  return(results)
}

parse_paf_asmgene <- function(paf_stats) {
  # Use lapply to read each file and convert it into a tibble
  results <- lapply(paf_stats, function(paf_stat) {
    haplotype <- str_extract(basename(paf_stat), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(paf_stat, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(paf_stat, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    
    # Read the file as a table and transform to long format
    read_tsv(paf_stat, col_names = c("col_type", "metric", "ref", "asm"), skip = 1) %>%
      select(-col_type) %>%
      pivot_longer(cols = c(ref, asm), names_to = "col", values_to = "value") %>%
      mutate(metric = str_c(col, ".", str_trim(metric)),  # Add prefix to metric
             value = case_when(         
               grepl("%$", value) ~ as.numeric(sub("%", "", value)) / 100,
               TRUE ~ suppressWarnings(as.numeric(value))
             ),
             haplotype = haplotype,
             asm_method = asm_method,
             asm_name = asm_name,
             source = "asmgene") %>%
      select(-col)
  })
  
  # Set names for the results list based on the input filenames
  names(results) <- paf_stats
  return(results)
}

parse_sex <- function(file) {
  results <- lapply(file, function(file) {
    haplotype <- str_extract(basename(file), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(file, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(file, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    
    read_tsv(file, col_names = c("metric", "value")) %>%
      mutate(metric = metric,
             value = as.numeric(value),
             haplotype = "both",
             asm_method = asm_method,
             asm_name = asm_name,
             source = "sex_determination")
  })
  names(results) <- file
  return(results)
}

parse_whatshap_stats <- function(file) {
  results <- lapply(file, function(file) {
    haplotype <- str_extract(basename(file), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(file, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(file, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    read_tsv(file, col_names = TRUE) %>%
      select(c(-`#sample`, -file_name)) %>%
      pivot_longer(cols = c(-chromosome), names_to = "metric", values_to = "value") %>%
      mutate(metric = metric,
             value = as.numeric(value),
             haplotype = "both",
             asm_method = asm_method,
             asm_name = asm_name,
             source = "whatshap_stats")
  })
  names(results) <- file
  return(results)
}

parse_whatshap_compare <- function(file) {
  results <- lapply(file, function(file) {
    haplotype <- str_extract(basename(file), "haplotype\\d")  # Extract haplotype from filename
    asm_method <- str_extract(file, "([^/]+)(?=/[^/]+/[^/]+$)")  # Extract assembly method from path
    asm_name <- str_extract(file, "([^/]+)(?=/[^/]+$)")  # Extract assembly name from path
    read_tsv(file, col_names = TRUE) %>%
      select(c(-`#sample`, -dataset_name0, -dataset_name1, -file_name0, -file_name1, -all_switchflips, -largestblock_switchflips)) %>%
      pivot_longer(cols = c(-chromosome), names_to = "metric", values_to = "value") %>%
      mutate(metric = metric,
             value = as.numeric(value),
             haplotype = "both",
             asm_method = asm_method,
             asm_name = asm_name,
             source = "whatshap_compare")
  })
  names(results) <- file
  return(results)
}

if (test_mode) {
  View(bind_rows(parse_paf_stat(input_paf_stat)))
  View(bind_rows(parse_paf_asmstat(input_paf_asmstat)))
  View(bind_rows(parse_paf_asmgene(input_paf_asmgene)))
  View(bind_rows(parse_sex(input_sex)))
  View(bind_rows(parse_whatshap_stats(input_whatshap_stats)))
  View(bind_rows(parse_whatshap_compare(input_whatshap_compare)))
}

bind_rows(
  bind_rows(parse_paf_stat(input_paf_stat)),
  bind_rows(parse_paf_asmstat(input_paf_asmstat)),
  bind_rows(parse_paf_asmgene(input_paf_asmgene)),
  bind_rows(parse_sex(input_sex)),
  bind_rows(parse_whatshap_stats(input_whatshap_stats)),
  bind_rows(parse_whatshap_compare(input_whatshap_compare))
) %>%
  write_tsv(output_qc_full)

