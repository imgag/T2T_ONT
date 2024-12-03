library(tidyverse)

plot_read_histograms <- function(length_files, qual_files, out_prefix) {
  # Read input files
  lengths <- list()
  for (file in length_files) {
    dt <- read_tsv(file,
      col_types = ("ddd")
    )
    filename <- basename(dirname(file))
    dt <- dt %>%
      add_column(filename) %>%
      set_names(c("lower_length", "upper_length", "count", "filename"))
    lengths[[filename]] <- dt
  }

  quals <- list()
  for (file in qual_files) {
    dt <- read_tsv(file,
      col_types = ("ddd")
    )
    filename <- basename(dirname(file))
    dt <- dt %>%
      add_column(filename) %>%
      set_names(c("lower_qual", "upper_qual", "count", "filename"))
    quals[[filename]] <- dt
  }

  l <- bind_rows(lengths)
  q <- bind_rows(quals)

  # Calculate N50 and total yield
  summary <- bind_cols(
    l %>%
      group_by(filename) %>%
      mutate(cumsum = cumsum(count * lower_length)) %>%
      mutate(N50_reached = cumsum > max(cumsum) / 2) %>%
      filter(N50_reached) %>%
      summarise(yield = max(cumsum),N50 = first(lower_length))
    ,
    l %>% 
      group_by(filename) %>%
      mutate(bases = count * lower_length) %>%
      summarise(
        longest_read = max(lower_length),
        bases_longer_1e5 = sum(bases[lower_length >= 1e5]),
        bases_longer_2e5 = sum(bases[lower_length >= 2e5]),
        reads_longer_1e6 = length(filename[lower_length >= 1e6]),
        ) %>%
      select(-filename)
  )
  
  write_tsv(summary, paste0("doc/tables/", out_prefix, ".summary.tsv"))

  # Plot length histogram bases_sequenced = count * length
  p1 <- ggplot(l) +
    stat_bin(aes(x = lower_length, weight = count, fill = filename),
      col = "grey20", position = "dodge"
    ) +
    theme_classic() +
    geom_vline(
      data = summary,
      aes(xintercept = N50, colour = filename), lty = 2
    )
  ggsave(paste0("doc/img/", out_prefix, ".read_length_hist.png"), plot = p1, width = 12, height = 6)

  # Plot Q histogram
  p2 <- ggplot(q) +
    stat_bin(aes(x = lower_qual, weight = count, fill = filename), col = "grey20", position = "dodge") +
    theme_classic()
  ggsave(paste0("doc/img/", out_prefix, ".read_qual_hist.png"), plot = p2, width = 12, height = 8)

  ## Combine both datasets into a single table
  long_l <- l %>%
    mutate(length = lower_length) %>%
    select(-c(upper_length, lower_length)) %>%
    pivot_longer(length, names_to = "measure")
  long_q <- q %>%
    mutate(qual = lower_qual) %>%
    select(-c(lower_qual, upper_qual)) %>%
    pivot_longer(qual, names_to = "measure")
  d <- bind_rows(long_l, long_q)
  
  # Plot combined histograms for qual and length
  p3 <- ggplot(d) +
    stat_bin(aes(x = value, weight = count, fill = filename),
      col = "grey20",
      position = "dodge"
    ) +
    theme_classic() +
    facet_wrap(~measure, scales = "free")
  ggsave(paste0("doc/img/", out_prefix, ".read_stats.png"), plot = p3, width = 12, height = 6)
}

# Compare UL published and sequenced
plot_read_histograms(
  length_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04399/run_04399.UL/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04400/run_04400.UL/length.hist"
  ),
  qual_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04399/run_04399.UL/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04400/run_04400.UL/quality.hist"
  ),
  out_prefix = "UL_comparison"
)

# Compare subsampling of UL published
plot_read_histograms(
  length_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.90x/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.70x/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.50x/length.hist"
  ),
  qual_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.90x/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.70x/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.50x/quality.hist"
  ),
  out_prefix = "UL_published_subsampling"
)

# Compare subsampling of HERRO
plot_read_histograms(
  length_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.120x/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.60x/length.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.35x/length.hist"
  ),
  qual_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.120x/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.60x/quality.hist",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.35x/quality.hist"
  ),
  out_prefix = "HQ_herro_published_subsampling"
)
