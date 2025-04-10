process_qc_table <- function(dt){

    # Calculate MMC
    dt_mmc <- dt %>%
    # First calculate MMC as before
    filter(metric %in% c("ref.dup_cnt", "asm.dup_cnt")) %>%
    pivot_wider(
        names_from = metric,
        values_from = value,
        values_fill = list(value = 0)
    ) %>%
    group_by(asm_name, haplotype, source, asm_method) %>%
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
    group_by(asm_name, haplotype, source, asm_method) %>%
    summarise(
        complete = sum(value[metric == "asm.full_sgl"]),
        total = sum(value[metric == "ref.full_sgl"]),
        value = (complete/total) * 100,
        .groups = 'drop'
    ) %>%
    select( -complete, -total) %>%
    mutate(metric = "genome_completeness")

    #head(dt_completeness)

    # add MMC and completeness to original table

    dt <- bind_rows(dt, dt_completeness, dt_mmc)

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


    # Remove chromosomes from whatshap stats
    # dt_reduced <- dt_reduced %>% 
    #     filter(!(str_detect(source, "whatshap_stats") & chromosome != "ALL")) %>%
    #     group_by(seq_type, metric, haplotype,source, asm_name) %>%
    #     reframe(
    #         value = case_when(
    #             str_detect(metric, "Rate|rate") ~ mean(value, na.rm = TRUE),
    #             metric == "Variants" ~ sum(value, na.rm = TRUE),
    #             TRUE ~ sum(value, na.rm = TRUE)
    #         ),
    #         .groups = 'drop'
    #     )

        return(dt_reduced)
}


