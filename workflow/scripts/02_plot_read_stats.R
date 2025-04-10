library(tidyverse)

### Plot density plots from readstats
plot_read_stats <- function(stats_files, out_prefix, save_plot = TRUE) {
  t <- list()
  for (f in stats_files) {
    t[[f]] <- read_tsv(f,
      col_select = c("sample_name", "read_length", "mean_quality")
      ) %>%
    slice_sample(n = 1e6)
  }

  dt <- bind_rows(t)
  "Finished data parsing"

  # Density plot
  p1 <- ggplot(dt, aes(x = read_length, y = mean_quality, col = sample_name)) +
    geom_density_2d() +
    ylim(c(0, 50)) +
    theme_classic()

  if (save_plot){
    plot_file = paste0("doc/img/", out_prefix, ".readstats_density.png")
    print("Writing plot to folder", plot_file)
    ggsave(plot_file, plot = p1)
  }
  return(p1)
}

