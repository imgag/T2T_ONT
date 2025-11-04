suppressPackageStartupMessages({
    library("optparse")   
    library("readr")
    library("dplyr")
    library("ggplot2")
    library("tidyr")})

# --------  input options  ------------------------------------
option_list <- list(
    make_option(c("--adj_insu"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--nonadj_insu"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--prefixes"), type = "character", 
    default = NULL, metavar = "character"),
    make_option(c("--plot"), type = "character", 
    default = NULL, metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$adj_insu) || is.null(opt$nonadj_insu)  || is.null(opt$plot) || is.null(opt$prefixes)) {
  stop("Input files must be specified")
}

# 1. --------------------------- plot all vs all insulation --------------

adj <- read_tsv(opt$adj_insu, show_col_types = F)[,c(1,2,6,8,10,12)] %>% 
    rename(`3x` = 3, `5x` = 4, `10x` = 5, `25x` = 6)

nonadj <- read_tsv(opt$nonadj_insu, show_col_types = F)[,c(1,2,6,8,10,12)] %>% 
    rename(`3x` = 3, `5x` = 4, `10x` = 5, `25x` = 6)

insulation <- merge(adj, nonadj, by = c("chrom","start"), suffixes = c("_adj","_nonadj")) %>%
    drop_na() %>%
    pivot_longer(
        cols = -c("chrom","start"),
        names_to = c("window","exp"),
        names_sep = "_",
        values_to = "insulation"
    ) %>%
    pivot_wider(
        names_from = exp,
        values_from = insulation
    ) %>%
    mutate(window = factor(window, levels = c("3x","5x","10x","25x")))

pcc <- insulation %>% group_by(window) %>%
    summarise(pcc = cor(adj, nonadj, use = "complete.obs", method="pearson"), .groups = "drop")

p1 <- ggplot(insulation, aes(x = adj, y = nonadj)) + 
    geom_point(alpha=0.5,size=0.4) +
    facet_wrap(~ window)+
    coord_cartesian(ylim=c(-4,3),xlim=c(-4,3))+
    labs(
        x="adj insulation",
        y="nonadj insulation"
    ) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    geom_text(
        data = pcc,
        aes(x=-3, y=2.5, 
        label = paste0("PCC=",signif(pcc,2))),
        size = 3, inherit.aes = F
    )

ggsave(filename = opt$plot, plot = p1, dpi = 300, width = 5, height = 5, units = "in")

#2 -------------------------- only plot bins that are boundaries in non adj -----

adj <- read_tsv(opt$adj_insu, show_col_types = F)[,c(1,2,6,8,10,12,18,19,20,21)] %>% 
rename(`3x` = 3,`5x` = 4,`10x` = 5,`25x` = 6, 
       `3x_boundary` = 7,`5x_boundary` = 8,`10x_boundary` = 9,`25x_boundary` = 10)

nonadj <- read_tsv(opt$nonadj_insu, show_col_types = F)[,c(1,2,6,8,10,12,18,19,20,21)] %>% 
  rename(`3x` = 3,`5x` = 4,`10x` = 5,`25x` = 6, 
         `3x_boundary` = 7,`5x_boundary` = 8,`10x_boundary` = 9,`25x_boundary` = 10)

insulation <- merge(adj, nonadj, by = c("chrom","start"), suffixes = c("_adj", "_nonadj"))

boundaries <- c("3x_boundary_nonadj", "5x_boundary_nonadj", 
                "10x_boundary_nonadj", "25x_boundary_nonadj")

for (boundary_col in boundaries) {  
  
  boundary_bins <- insulation %>% filter(.data[[boundary_col]] == TRUE) %>%
    select(-contains("boundary")) %>%
    drop_na() %>%
    pivot_longer(
      cols = -c(chrom, start),
      names_to = c("window", "exp"),
      names_sep = "_",
      values_to = "insulation"
    ) %>%
    pivot_wider(
      names_from = exp,
      values_from = insulation
    ) %>% mutate(window = factor(window,levels=c("3x","5x","10x","25x")))
  
  pcc <- boundary_bins %>%
    group_by(window) %>%
    summarise(pcc = cor(adj, nonadj,use = "complete.obs", method = "pearson"), .groups = "drop")
  
  p2 <- ggplot(boundary_bins, aes(x = adj, y = nonadj)) +
    geom_point(alpha = 0.6, size = 0.4) +
    facet_wrap(~ window) +
    coord_cartesian(ylim = c(-4,3),xlim=c(-4,3)) +
    labs(
      x = "insulation score adj",
      y = "insulation score nonadj",
      title = paste0("is boundary in: " ,boundary_col)) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    geom_text(
      data = pcc,
      aes(x = -3, y = 2.5, label = paste0("PCC = ", signif(pcc, 2))), size = 3,
      inherit.aes = FALSE
    )
     
  ggsave(filename = paste0(opt$prefixes, boundary_col, ".png"), 
  plot = p2, dpi = 300, width = 5, height = 5, units = "in")

}

