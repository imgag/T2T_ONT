library(tidyverse)

### Plot density plots from readstats
plot_read_stats <- function(stats_files, out_prefix) {
  t <- list()
  for (f in fl) {
    t[[f]] <- read_tsv(f,
      col_select = c("sample_name", "read_length", "mean_quality")
    ) %>%
      slice_sample(n = 1e5)
  }

  dt <- bind_rows(t)

  # Density plot
  p1 <- ggplot(dt, aes(x = read_length, y = mean_quality, col = sample_name)) +
    geom_density_2d() +
    ylim(c(0, 50)) +
    theme_classic()

  ggsave(paste0("doc/img/", out_prefix, ".readstats_density.png", plot = p1))
}


# Plot read length and qual densities of the two UL flowcells, along with published UL and PoreC and Duplex for context
plot_read_stats(
  stats_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_duplex/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.POREC/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04399/run_04399.UL/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04400/run_04400.UL/read_stats.txt"
  ),
  out_prefix = "UL_full_context"
)

# Plot read length and qual densities for downsamples UL flowcells (published)
plot_read_stats(
  stats_files = c(
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL50x/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL70x/read_stats.txt",
    "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL90x/read_stats.txt",
    ),
  out_prefix = "UL_subsampling"
)