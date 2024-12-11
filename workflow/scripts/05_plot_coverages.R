library(tidyverse)

plot_boxplots <- function(cov_files, out_prefix) {

  l <- list()
  for (f in cov_files) {
    d <- read_tsv(f,
      col_names = c("chr", "cov", "percentile")
    ) %>%
      filter(chr == "total") %>%
      add_column(sample = basename(dirname(f)))  %>%
      separate(sample, 
               into = c("dataset", "method", "coverage"),
               sep = "\\.",
               fill = "right") %>%
      add_column(sample = basename(dirname(f)))
    l[[f]] <- d
  }

  dt <- bind_rows(l)

  dt <- dt %>% 
    mutate(coverage = ifelse(is.na(coverage), "all", coverage)) %>% 
    arrange(dataset, method, coverage)

  dt <- dt %>%
    group_by(sample, dataset, method, coverage) %>%
    summarise(
      Q1 = max(cov[percentile >= 0.75]),
      median = max(cov[percentile >= 0.5]),
      Q3 = max(cov[percentile >= 0.25]),
      lower_whisker = max(Q1 - (1.5 * (Q3 - Q1)), 0),
      upper_whisker = max(Q3 + (1.5 * (Q3 - Q1)), 0),
    )
  
  p <- ggplot(dt, aes(x = paste(method, coverage), fill = method)) +
    geom_boxplot(aes(
      ymin = lower_whisker,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = upper_whisker
    ), stat = "identity") +
    #facet_wrap(~dataset, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Mapping coverage") +
    theme(axis.title.x = element_blank()) 

  ggsave(paste0("doc/img/", out_prefix, ".cov_boxplot.png"), plot = p, width = 16, height = 6)
}

plot_boxplots(
  cov_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_duplex/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_duplex.15x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_duplex.20x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_duplex.25x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_herro/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_herro.35x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_herro.60x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_herro.120x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.UL/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.UL.50x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.UL.70x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.UL.90x/cov.mosdepth.global.dist.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_combined.15x_20x/cov.mosdepth.global.dist.txt"

  ), 
  out_prefix = "subsampling_coverages"
)
