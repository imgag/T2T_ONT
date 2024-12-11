library(tidyverse)

### Comparison Mapping stats of selected datasets in a single plot 

f1 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_duplex/bamstats.txt"
f2 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.HQ_herro/bamstats.txt"
f3 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/published/published.UL/bamstats.txt"
f4 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/run_04399/run_04399.UL/bamstats.txt"
f5 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/bamstats/run_04400/run_04400.UL/bamstats.txt"

fl <- list(f1, f2, f3, f4, f5)
t <- list()
for (f in fl){
  t[[f]] <- read_tsv(f)
}
dt <- bind_rows(t)

dtp <- dt %>%
    mutate(sample_name = as_factor(sample_name)) %>%
    group_by(sample_name) %>%
    slice_sample(n=1e5) %>%
    select(-c(runid, start_time)) %>%
    drop_na()
    
dtp %>% head() %>% view()

p2 <- ggplot(dtp %>% filter(sample_name != "published/published.HQ_herro"),
        aes(x=read_length, y=mean_quality, col = sample_name)) +
    geom_density_2d(contour_var = "count") +
    theme_classic()

ggsave("doc/img/mapstats_len_quality.png", plot = p2)

p3 <- ggplot(dtp, aes(x=sample_name, y=mean_quality, col = sample_name) ) +
  geom_boxplot() +
  theme_classic()
ggsave("doc/img/mapstats_qual_boxplot.png", plot = p3)
p3

p4 <- ggplot(dtp, aes(x=read_length, col = sample_name) ) +
  geom_density() +
  theme_classic() +
  scale_x_log10()
ggsave("doc/img/mapstats_read_length_density.png", plot = p4)
p4
