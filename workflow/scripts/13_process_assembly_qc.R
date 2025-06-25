selected_metrics = c(
    "length",
    "ref_covered",
    "asm_covered",
    "inter-chromosomal_misjoins",
    "intra-chromosomal_gaps",
    "candidate_inversions_in_the_middle",
    "candidate_inversions_at_contig_ends",
    "variants",
    "blocks",
    "qv",
    "n_T2T",
    "total_n_count",
    "total_gaps",
    "n_contigs",
    "n_contigs_over_10mb",
    "genome_completeness",
    "MMC",
    "overall_switch_rate"
    )

process_qc_table <- function(dt, selection = selected_metrics) {

    # Define the optional grouping columns
    optional_group_cols <- c("n_UL", "n_DX", "sample")
    
    # Get the base grouping columns that are always present
    base_group_cols <- c("asm_name", "haplotype", "source", "asm_method")
    
    # Add optional columns if they exist in the data
    group_cols <- c(base_group_cols, optional_group_cols[optional_group_cols %in% names(dt)])
    
    # Calculate MMC
    dt_mmc <- dt %>%
    # First calculate MMC as before
    filter(metric %in% c("ref.dup_cnt", "asm.dup_cnt")) %>%
    pivot_wider(
        names_from = metric,
        values_from = value,
        values_fill = list(value = 0)
        ) %>%
        group_by(!!!syms(group_cols)) %>%
    summarise(
        ref_mc = sum(ref.dup_cnt),
        asm_mc = sum(asm.dup_cnt),
        value = 1 - asm_mc/ref_mc,
        .groups = 'drop'
    ) %>% 
    select( -ref_mc, -asm_mc) %>%
    mutate(metric = "MMC")

    #head(dt_mmc)

    # Calculate genome completeness
    dt_completeness <- dt %>%
    filter(source == "asmgene") %>%
    group_by(!!!syms(group_cols)) %>%
    summarise(
        complete = sum(
            value[metric == "asm.full_sgl"],
            value[metric == "asm.full_dpl"]
            ),
        total = sum(value[metric == "ref.full_sgl"]),
        value = (complete/total) * 100,
        .groups = 'drop'
    ) %>%
    select( -complete, -total) %>%
    mutate(metric = "genome_completeness")

    dt_multiple_genes <- dt %>%
    filter(source == "asmgene") %>%
    group_by(!!!syms(group_cols)) %>%
    summarise(
        complete = sum(
            value[metric == "asm.full_sgl"],
            value[metric == "asm.full_dpl"]
            ),
        total = sum(value[metric == "ref.full_sgl"]),
        value = (complete/total) * 100,
        .groups = 'drop'
    ) %>%
    select( -complete, -total) %>%
    mutate(metric = "genome_completeness")
    #head(dt_completeness)

    # Aggregate whatshap_compare metrics across chromosomes
    dt_whatshap_agg <- dt %>%
    filter(source == "whatshap_compare") %>%
    group_by(!!!syms(group_cols)) %>%
    summarise(
        # Calculate total covered variants
        total_covered_variants = sum(value[metric == "covered_variants"], na.rm = TRUE),
        
        # Calculate weighted switch rate
        weighted_switch_rate = sum(
            value[metric == "covered_variants"] * value[metric == "all_switch_rate"] / 100, 
            na.rm = TRUE
        ),
        
        # Calculate weighted switch flip rate
        weighted_switch_flip_rate = sum(
            value[metric == "covered_variants"] * value[metric == "all_switchflip_rate"] / 100, 
            na.rm = TRUE
        ),
        .groups = 'drop'
    ) %>%
    # Calculate the overall rates based on weightings
    mutate(
        overall_switch_rate = weighted_switch_rate / total_covered_variants * 100,
        overall_switch_flip_rate = weighted_switch_flip_rate / total_covered_variants * 100
    ) %>%
    select(-weighted_switch_rate, -weighted_switch_flip_rate) %>%
    # Reshape to long format
    pivot_longer(
        cols = c(total_covered_variants, overall_switch_rate, overall_switch_flip_rate),
        names_to = "metric",
        values_to = "value"
    ) %>%
    mutate(
        source = "whatshap_compare_aggregated",
        chromosome = "ALL"  # Marking as aggregated across all chromosomes
    )

    # add MMC, completeness, and whatshap_agg to original table
    dt <- bind_rows(dt, dt_completeness, dt_mmc, dt_whatshap_agg)

    key_metrics = c(
        "")
    # Remove  values from dataset
    dt_reduced <- dt %>%
        filter(metric %in% selection)

    return(dt_reduced)
}



