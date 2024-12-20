library(tidyverse)


files <- list()

parse_summary_summaries <- function (files, out_f){

    dt <- bind_rows(
        lapply(files, function(f) {
            s <- str_split(basename(f), "[.]", simplify = T)
            run <- s[1]
            type <- s[2]

            d <- read_tsv(f) %>%
                add_column("run" = run) %>%
                add_column("type" = type) %>%
                add_column("filepath" = f) %>%
                relocate()

            return(d)
        })
    )

    write_tsv(dt, out_f)
}

UL_files <- list.files(
    path="data/bamstats/basecalled/SUP",
    pattern="*.sequencing_summary.processed.txt",
    recursive = T,
    full.names = T)

parse_summary_summaries(UL_files, "doc/tables/UL_sequencing_stats.tsv")
