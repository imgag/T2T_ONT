#!/usr/bin/env Rscript

# Script to aggregate values from flowcells_qc.tsv where source == whatshap_compare
# Sum up the covered variants and recalculate rates

library(tidyverse)

# Read the data
data <- read_tsv("doc/tables/flowcells_qc.tsv")

# Filter for whatshap_compare rows
whatshap_data <- data %>%
  filter(source == "whatshap_compare")

# Group by key identifying columns
aggregated_data <- whatshap_data %>%
  group_by(n_UL, n_DX, sample, asm_name, haplotype, asm_method) %>%
  summarize(
    # Calculate total covered variants
    total_covered_variants = sum(value[metric == "Covered Variants"], na.rm = TRUE),
    
    # Calculate weighted switch rate
    weighted_switch_rate = sum(
      value[metric == "Covered Variants"] * value[metric == "Switch Rate (%)"] / 100, 
      na.rm = TRUE
    ),
    
    # Calculate weighted switch flip rate
    weighted_switch_flip_rate = sum(
      value[metric == "Covered Variants"] * value[metric == "Switch Flip Rate (%)"] / 100, 
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  # Calculate the actual rates based on weightings
  mutate(
    overall_switch_rate = weighted_switch_rate / total_covered_variants * 100,
    overall_switch_flip_rate = weighted_switch_flip_rate / total_covered_variants * 100
  ) %>%
  # Clean up intermediate columns
  select(-weighted_switch_rate, -weighted_switch_flip_rate)

# Create tidy format for the results
tidy_results <- aggregated_data %>%
  pivot_longer(
    cols = c(total_covered_variants, overall_switch_rate, overall_switch_flip_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_when(
      metric == "total_covered_variants" ~ "Covered Variants (Total)",
      metric == "overall_switch_rate" ~ "Overall Switch Rate (%)",
      metric == "overall_switch_flip_rate" ~ "Overall Switch Flip Rate (%)",
      TRUE ~ metric
    ),
    source = "whatshap_compare_aggregated",
    chromosome = "ALL"  # Marking as aggregated across all chromosomes
  )

# Print results
print(tidy_results)

# Write results to file
write_tsv(tidy_results, "doc/tables/flowcells_qc_whatshap_aggregated.tsv")

# Also print a nice summary table
summary_table <- tidy_results %>%
  select(n_UL, n_DX, sample, asm_name, metric, value) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

print(summary_table) 