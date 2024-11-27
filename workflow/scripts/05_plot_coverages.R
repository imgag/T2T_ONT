library(tidyverse)

plot_boxplots <- function(cov_files, out_prefix) {
  l <- list()
  for (f in files) {
    d <- read_tsv(f,
      col_names = c("chr", "cov", "percentile")
    ) %>%
      filter(chr == "total") %>%
      add_column(sample = basename(dirname(f)))
    l[[f]] <- d
  }

  dt <- bind_rows(l)

  dt <- dt %>%
    group_by(sample) %>%
    summarise(
      Q1 = max(cov[percentile >= 0.75]),
      median = max(cov[percentile >= 0.5]),
      Q3 = max(cov[percentile >= 0.25]),
      lower_whisker = max(Q1 - (1.5 * (Q3 - Q1)), 0),
      upper_whisker = max(Q3 + (1.5 * (Q3 - Q1)), 0),
    )

  p <- ggplot(dt, aes(x = sample, fill = sample)) +
    geom_boxplot(aes(
      ymin = lower_whisker,
      lower = Q1,
      middle = median,
      upper = Q3,
      ymax = upper_whisker
    ), stat = "identity") +
    theme_classic()

  ggsave(paste0("doc/img/", prefix, ".cov_boxplot.png"), plot = p)
}

# Print coverages for subsampled Herro reads
# todo calculate all coverages!

plot_boxplots(
  cov_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/run_04399/run_04399.UL/cov.mosdepth.global.dist.txt"
  )
)
