suppressPackageStartupMessages({
    library("optparse")   
    library("readr")
    library("dplyr")
    library("ggplot2")
    library("tidyr")})

# --------  input options  ------------------------------------
option_list <- list(
    make_option(c("--frags"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--plot"), type = "character", 
    default = NULL, metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$frags)) {
  stop("Input files must be specified")
}

# ----- Plot all vs all compartment scores (eig 1, pos = A, neg = B)
frags <- read.table(opt$frags, col.names = c("reads","frags"))
cutoff = 21
summary_row <- data.frame(reads = sum(frags$reads[frags$frags >= cutoff]), frags = paste(cutoff, "+"))

frags <- rbind(frags[frags$frags < cutoff, ], summary_row)
frags$frags <- factor(frags$frags, levels = c(1:20,"21 +"))
#total_reads <- sum(frags$reads)
#total_frags <- sum(frags$frags*frags$reads)


p <- ggplot(frags, aes(x = as.factor(frags), y=reads)) +
    geom_bar(stat = "identity") +
    labs(
        x="number of contacts from 1 read",
        y="read count"
    ) +
    theme_bw() +
    theme(aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) 

ggsave(filename = opt$plot, plot = p, dpi = 300)