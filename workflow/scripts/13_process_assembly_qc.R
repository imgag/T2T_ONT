default_selected_metrics = c(
    "length",
    "n_contigs_over_10mb",
    "n_contigs",
    "ref_covered",
    "asm_covered",
    "genome_completeness",
    "MMC",
    "MSC",
    "inter-chromosomal_misjoins",
    "intra-chromosomal_gaps",
    "candidate_inversions_in_the_middle",
    "candidate_inversions_at_contig_ends",
    "n_t2t_contig",
    "n_t2t_scaffold",
    "total_gaps",
    "total_n_count",
    "total_gaps",
    "total_covered_variants",
    "overall_switch_rate",
    "overall_switch_flip_rate",
    "total_blocks",
    "total_covered_variants",
    "total_switches",
    "total_het_variants",
    "overall_switch_rate",
    "overall_switchflip_rate",
    "overall_hamming_rate",
    "overall_diff_genotypes_rate"
    )

process_qc_table <- function(dt, selection = default_selected_metrics) {

    # Define the optional grouping columns
    optional_group_cols <- c("n_UL", "n_DX", "sample")
    
    # Get the base grouping columns that are always present
    base_group_cols <- c("asm_name", "haplotype", "source", "asm_method")
    
    # Add optional columns if they exist in the data
    group_cols <- c(base_group_cols, optional_group_cols[optional_group_cols %in% names(dt)])
    
    # Calculate MMC (Missing Multi Copy)
    dt_mmc <- dt %>%
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

    # Calculate MSC (Missing Single Copy)
    dt_msc <- dt %>%
    filter(source == "asmgene") %>%
    group_by(!!!syms(group_cols)) %>%
    summarise(
        ref_single = value[metric == "ref.full_sgl"],
        asm_single = value[metric == "asm.full_sgl"],
        value = (ref_single - asm_single) / ref_single,
        .groups = 'drop'
    ) %>%
    select(-ref_single, -asm_single) %>%
    mutate(metric = "MSC")

    
    # Calculate genome completeness
    dt_completeness <- dt %>%
    filter(source == "asmgene") %>%
    group_by(!!!syms(group_cols)) %>%
    summarise(
        complete = sum(
            value[metric == "asm.full_sgl"],
            value[metric == "asm.full_dup"]
            ),
        total = sum(value[metric == "ref.full_sgl"]),
        value = (complete/total),
        .groups = 'drop'
    ) %>%
    select( -complete, -total) %>%
    mutate(metric = "genome_completeness")

    # Aggregate whatshap_compare metrics across chromosomes
    # First pivot to wide format per chromosome for easier calculations
    dt_whatshap_wide <- dt %>%
        filter(source == "whatshap_compare") %>%
        pivot_wider(
            names_from = metric,
            values_from = value,
            values_fill = list(value = 0)
        )
    
    # Aggregate metrics across chromosomes
    dt_whatshap_agg <- dt_whatshap_wide %>%
        group_by(!!!syms(group_cols)) %>%
        summarise(
            # Count metrics (sum)
            total_blocks = sum(intersection_blocks, na.rm = TRUE),
            total_covered_variants = sum(covered_variants, na.rm = TRUE),
            total_switches = sum(all_switches, na.rm = TRUE),
            total_het_variants = sum(het_variants0, na.rm = TRUE),
            total_blockwise_hamming = sum(blockwise_hamming, na.rm = TRUE),
            total_blockwise_diff_genotypes = sum(blockwise_diff_genotypes, na.rm = TRUE),
            
            # Weighted rates (weighted by covered_variants)
            weighted_switch_rate = sum(covered_variants * all_switch_rate, na.rm = TRUE),
            weighted_switchflip_rate = sum(covered_variants * all_switchflip_rate, na.rm = TRUE),
            weighted_hamming_rate = sum(covered_variants * blockwise_hamming_rate, na.rm = TRUE),
            weighted_diff_genotypes_rate = sum(covered_variants * blockwise_diff_genotypes_rate, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        # Calculate the overall rates based on weightings
        mutate(
            overall_switch_rate = weighted_switch_rate / total_covered_variants,
            overall_switchflip_rate = weighted_switchflip_rate / total_covered_variants,
            overall_hamming_rate = weighted_hamming_rate / total_covered_variants,
            overall_diff_genotypes_rate = weighted_diff_genotypes_rate / total_covered_variants
        ) %>%
        select(-weighted_switch_rate, -weighted_switchflip_rate, 
               -weighted_hamming_rate, -weighted_diff_genotypes_rate,
               -total_blockwise_hamming, -total_blockwise_diff_genotypes) %>%
        # Reshape to long format
        pivot_longer(
            cols = c(total_blocks, total_covered_variants, total_switches, total_het_variants,
                     overall_switch_rate, overall_switchflip_rate, 
                     overall_hamming_rate, overall_diff_genotypes_rate),
            names_to = "metric",
            values_to = "value"
        ) %>%
        mutate(
            source = "whatshap_compare_aggregated",
            chromosome = "ALL"  # Marking as aggregated across all chromosomes
        )

    # add MMC, completeness, and whatshap_agg to original table
    dt <- bind_rows(dt, dt_completeness, dt_mmc, dt_msc, dt_whatshap_agg)

    key_metrics = c(
        "")
    # Remove  values from dataset
    dt_reduced <- dt %>%
        filter(metric %in% selection)

    return(dt_reduced)
}



