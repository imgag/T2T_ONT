suppressPackageStartupMessages({
    library("optparse")   
    library("readr")
    library("dplyr")
    library("ggplot2")
    library("tidyr")})

# --------  input options  ------------------------------------
option_list <- list(
    make_option(c("--adj_eigs"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--nonadj_eigs"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--plot"), type = "character", 
    default = NULL, metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$adj_eigs) || is.null(opt$nonadj_eigs)  || is.null(opt$plot)) {
  stop("Input files must be specified")
}

# ----- Plot all vs all compartment scores (eig 1, pos = A, neg = B)
adj <- read_tsv(opt$adj_eigs, show_col_types = F)[,c("chrom","start","E1")] 
nonadj <- read_tsv(opt$nonadj_eigs, show_col_types = F)[,c("chrom","start","E1")] 

eigs <- merge(adj, nonadj, by = c("chrom","start"), suffixes = c("_adj","_nonadj")) %>%
    drop_na()

pcc <- eigs %>%
    group_by(chrom) %>%
    summarise(pcc = cor(E1_adj,E1_nonadj,use="complete.obs",method="pearson"), .groups = "drop")

p <- ggplot(eigs, aes(x = E1_adj, y=E1_nonadj)) +
    geom_point(alpha=0.5,size=0.3) +
    facet_wrap(~ chrom)+
    coord_cartesian(ylim=c(-2,2),xlim=c(-2,2)) + 
    labs(
        x="adj eig 1",
        y="nonadj eig 1"
    ) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    geom_text(
        data = pcc,
        aes(x=-0.35, y=1.3, 
        label = paste0("PCC=",signif(pcc,2))),
        size = 3, inherit.aes = F
    )

ggsave(filename = opt$plot, plot = p, dpi = 300)