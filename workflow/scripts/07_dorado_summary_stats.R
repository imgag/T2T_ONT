library(tidyverse)

log_file <- snakemake@log[[1]]
sink(log_file, append = TRUE, type = "output")
sink(log_file, append = TRUE, type = "message")

f <- snakemake@input[[1]]

dt <- read_tsv(f, col_names = T) %>%
    select(c("mean_qscore_template", "sequence_length_template", "start_time", "duration", "mux")) %>%
    set_names(c("qscore", "length", "start_time", "duration", "mux"))
summary <- bind_cols(
    dt %>% 
        summarise(
            yield = sum(length),
            longest_read = max(length),
            bases_longer_1e5 = sum(length[length >= 1e5]),
            bases_longer_2e5 = sum(length[length >= 2e5]),
            reads_longer_1e6 = length(length[length >= 1e6]),
        ),
    dt %>%
        mutate(cumsum = cumsum(length)) %>%
        mutate(N50_reached = cumsum > max(cumsum)/2) %>%
        filter(N50_reached) %>%
        summarise(N50 = first(cumsum)) 
)

write_tsv(summary, snakemake@output[[1]])