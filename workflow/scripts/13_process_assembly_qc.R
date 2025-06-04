process_qc_table <- function(dt){

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
        metric = case_when(
            metric == "total_covered_variants" ~ "Covered Variants (Total)",
            metric == "overall_switch_rate" ~ "Overall Switch Rate (%)",
            metric == "overall_switch_flip_rate" ~ "Overall Switch Flip Rate (%)",
            TRUE ~ metric
        ),
        source = "whatshap_compare_aggregated",
        chromosome = "ALL"  # Marking as aggregated across all chromosomes
    )

    # add MMC, completeness, and whatshap_agg to original table
    dt <- bind_rows(dt, dt_completeness, dt_mmc, dt_whatshap_agg)

    # Remove  values from dataset
    dt_reduced <- dt %>%
        filter(!str_detect(metric, "Number of")) %>%
        filter(!str_detect(metric, "asm.")) %>% filter(!str_detect(metric, "ref.")) %>%
        filter(!(source == "merqury" & str_detect(metric, "_kmers"))) %>%
        filter(!(source == "whatshap_compare" & str_detect(metric, "block"))) %>%
        filter(!(source == "whatshap_stats" & str_detect(metric, "_"))) %>%
        filter(!(source == "whatshap_stats" & str_detect(metric, "_"))) %>%
        filter(!(source == "t2t_motif")) %>%
        filter(!(source == "whatshap_compare" & metric %in% c("het_variants0", "only_snvs", 'all_assessed_pairs', 'all_switches'))) %>% 
        filter(!(source == "whatshap_stats" &  metric %in% c("covered_variants", "only_snvs"))) %>%
        filter(!(source == "asmstat" & str_detect(metric, "bp"))) %>%
        filter(!(source == "asmstat" & str_detect(metric, "l_cov")))

        # Rename metrics to better names used in plotting
    dt_reduced <- dt_reduced %>%
        mutate(
            metric = case_when(
                metric == "NG50" ~ "NG50",
                metric == "n_T2T" ~ "Number of T2T chromosomes",
                metric == "qv" ~ "Quality Value",
                metric == "error_rate" ~ "Error Rate",
                metric == "all_switch_rate" ~ "Switch Rate (%)",
                metric == "#breaks" ~ "Number of Breaks in Assembly",
                metric == "MMC" ~ "Missing Multi-Copy Genes (%)",
                metric == "Length" ~ "Assembly Length",
                metric == "Rcov" ~ "% of Ref covered by Assembly",
                metric == "Rdup" ~ "% of Ref duplicated in Assembly",
                metric == "Qcov" ~ "% of Assembly covered by Ref",
                metric == "NG75" ~ "NG75",
                metric == "NGA50" ~ "NGA50",
                metric == "AUNGA" ~ "Average Ungapped Alignment length",
                metric == "variants" ~ "Variants",
                metric == "phased" ~ "Variants Phased",
                metric == "unphased" ~ "Variants Unphased",
                metric == "blocks" ~ "Phase Blocks",
                metric == "covered_variants" ~ "Covered Variants",
                metric == "all_switchflip_rate" ~ "Switch Flip Rate (%)",
                metric == "n_y_chrom" ~ "Number of Y Chromosomes",
                metric == "genome_completeness" ~ "Genome Completeness",
                TRUE ~ metric
                
            )
        )


    # Keep only ALL chromosome rows from whatshap_stats, remove per-chromosome stats
    dt_reduced <- dt_reduced %>%
        filter(!(source == "whatshap_stats" & chromosome != "ALL")) %>%
        filter(!(source == "whatshap_compare")) %>%
        filter(!(source == "whatshap_compare_aggregated" & chromosome != "ALL"))

    return(dt_reduced)
}



