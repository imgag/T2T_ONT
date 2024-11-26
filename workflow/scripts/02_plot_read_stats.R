library(tidyverse)

### Plot density plots from readstats

fl <- c(
  "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_duplex/read_stats.txt",
  "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL/read_stats.txt",
  "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.POREC/read_stats.txt",
  "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04399/run_04399.UL/read_stats.txt",
  "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/run_04400/run_04400.UL/read_stats.txt"
)

t <- list()
for (f in fl){
  t[[f]] <- read_tsv(f,
  col_select = c("sample_name", "read_length", "mean_quality")) %>%
  slice_sample(n=1e5)
}

dt <- bind_rows(t)

# Density plot
p1 <- ggplot(dt, aes(x=read_length, y=mean_quality, col = sample_name) ) +
  geom_density_2d() +
  ylim(c(0,50)) +
  theme_classic()

ggsave("doc/img/readstats_len_qual.png", plot = p1)

# Calculate bases sequenced

dt_bin <- dt %>%
  mutate(bin_quality = cut_interval(mean_quality, n = 10)) %>%
  mutate(bin_length = cut_interval(read_length, n = 100)) %>%
  group_by(bin_length) %>%
  mutate(bases_sequenced = sum(read_length)) %>%
  select(sample_name, bin_quality, bin_length, bases_sequenced) %>%
  unique()

p2 <- ggplot(dt_bin, aes(x = bin_length, y = bases_sequenced, fill = bin_quality)) +
  geom_bar()
