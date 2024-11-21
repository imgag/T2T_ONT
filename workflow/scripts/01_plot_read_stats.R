library(tidyverse)


length_files <- list.files("assembly/input_qc/published", recursive = TRUE, pattern = "length.hist", full.names=TRUE)
qual_files <- list.files("assembly/input_qc/published", recursive = TRUE, pattern ="quality.hist", full.names=TRUE)


lengths <- list()
for (file in length_files) {
  dt <- read_tsv(file)
  filename <- basename(dirname(file))
  dt <- dt %>%
    add_column(filename) %>%
    set_names(c("lower_length", "upper_length", "count", "filename"))
  lengths[[filename]] <- dt
}

quals <- list()
for (file in qual_files) {
  dt <- read_tsv(file)
  filename <- basename(dirname(file))
  dt <- dt %>% 
    add_column(filename) %>%
    set_names(c("lower_qual", "upper_qual", "count", "filename"))
  quals[[filename]] <- dt
}

all_lengths <- bind_rows(lengths)
all_lengths <- bind_cols(quals)

### Comparison Read stats of selected datasets in a single plot 

f1 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_duplex/read_stats.txt"
f2 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.UL.70x/read_stats.txt"
f3 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.POREC/read_stats.txt"
#f4 <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input_qc/published/published.HQ_herro.70x/read_stats.txt"
# run1
# run2

fl <- list(f1, f2, f3, f4)
t <- list()
for (f in fl){
  t[[f]] <- read_tsv(f,
  col_select = c("sample_name", "read_length", "mean_quality")) %>%
  slice_sample(n=1e5)
}

dt <- bind_rows(t)

p2 <- ggplot(dt, aes(x=read_length, y=mean_quality, col = sample_name) ) +
  geom_density_2d() +
  ylim(c(0,50)) +
  
  theme_classic()
p2

ggsave("doc/img/readstats_len_qual.png", plot = p2)


dt <- full_join(bind_rows(lengths), bind_rows(quals))
