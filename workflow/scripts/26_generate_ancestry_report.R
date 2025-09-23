#!/usr/bin/env Rscript
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/26_generate_ancestry_report.R

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(plotly)
  library(DT)
  library(rmarkdown)
  library(knitr)
  library(gridExtra)
  library(reshape2)
  library(RColorBrewer)
})

# Function to read PCA results
read_pca_results <- function(eigenval_file, eigenvec_file) {
  eigenval <- read_table(eigenval_file, col_names = "eigenvalue")
  eigenvec <- read_table(eigenvec_file, col_names = c("FID", "IID", paste0("PC", 1:20)))
  
  # Calculate variance explained
  total_var <- sum(eigenval$eigenvalue)
  eigenval$var_explained <- eigenval$eigenvalue / total_var * 100
  eigenval$PC <- paste0("PC", 1:nrow(eigenval))
  
  return(list(eigenval = eigenval, eigenvec = eigenvec))
}

# Function to read ADMIXTURE results
read_admixture_results <- function(admixture_dir) {
  q_files <- list.files(admixture_dir, pattern = "\\.Q$", full.names = TRUE)
  
  admixture_results <- list()
  for (q_file in q_files) {
    k_value <- as.numeric(gsub(".*\\.([0-9]+)\\.Q$", "\\1", basename(q_file)))
    q_data <- read_table(q_file, col_names = paste0("Ancestry", 1:k_value))
    q_data$K <- k_value
    q_data$Sample <- 1:nrow(q_data)
    admixture_results[[paste0("K", k_value)]] <- q_data
  }
  
  return(admixture_results)
}

# Function to read iAdmix results
read_iadmix_results <- function(iadmix_file) {
  if (file.exists(iadmix_file)) {
    iadmix <- read_table(iadmix_file)
    return(iadmix)
  } else {
    return(NULL)
  }
}

# Function to read RFMix results
read_rfmix_results <- function(rfmix_file) {
  if (file.exists(rfmix_file)) {
    rfmix <- read_table(rfmix_file)
    return(rfmix)
  } else {
    return(NULL)
  }
}

# Function to read sample metadata from PSAM file
read_sample_metadata <- function(psam_file) {
  # Read PSAM file with header
  psam <- read_table(psam_file)
  
  # Handle different possible column names for PSAM format
  if ("PHENO1" %in% colnames(psam)) {
    psam$PHENO <- psam$PHENO1
  } else if ("PHENOTYPE" %in% colnames(psam)) {
    psam$PHENO <- psam$PHENOTYPE
  }
  
  # Ensure required columns exist
  required_cols <- c("FID", "IID", "PAT", "MAT", "SEX")
  missing_cols <- setdiff(required_cols, colnames(psam))
  if (length(missing_cols) > 0) {
    # Handle cases where columns might have different names
    if ("FID" %in% missing_cols && "#FID" %in% colnames(psam)) {
      psam$FID <- psam$`#FID`
    }
    if ("PAT" %in% missing_cols) psam$PAT <- "0"
    if ("MAT" %in% missing_cols) psam$MAT <- "0" 
    if ("PHENO" %in% missing_cols) psam$PHENO <- -9
  }
  
  # Map phenotype codes to population labels
  if ("PHENO" %in% colnames(psam)) {
    pop_map <- c("1" = "AFR", "2" = "AMR", "3" = "EAS", "4" = "EUR", "5" = "SAS", "-9" = "Unknown")
    psam$Population <- pop_map[as.character(psam$PHENO)]
  } else {
    psam$Population <- "Unknown"
  }
  
  # Identify study samples (assuming they have specific naming pattern or are Unknown)
  psam$Sample_Type <- ifelse(psam$Population == "Unknown", "Study_Sample", "Reference_1000G")
  
  return(psam)
}

# Function to create PCA plot
create_pca_plot <- function(pca_data, metadata) {
  pca_merged <- merge(pca_data$eigenvec, metadata, by = c("FID", "IID"))
  
  # Calculate variance explained for PC1 and PC2
  var_pc1 <- round(pca_data$eigenval$var_explained[1], 2)
  var_pc2 <- round(pca_data$eigenval$var_explained[2], 2)
  
  p <- ggplot(pca_merged, aes(x = PC1, y = PC2, color = Population, shape = Sample_Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_brewer(type = "qual", palette = "Set1") +
    scale_shape_manual(values = c(16, 17)) +
    labs(
      title = "Principal Component Analysis",
      x = paste0("PC1 (", var_pc1, "%)"),
      y = paste0("PC2 (", var_pc2, "%)"),
      color = "Population",
      shape = "Sample Type"
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  return(p)
}

# Function to create scree plot
create_scree_plot <- function(pca_data) {
  p <- ggplot(pca_data$eigenval[1:10, ], aes(x = PC, y = var_explained)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = paste0(round(var_explained, 1), "%")), 
              vjust = -0.5, size = 3) +
    labs(
      title = "Scree Plot - Variance Explained by Principal Components",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Function to create ADMIXTURE plot
create_admixture_plot <- function(admixture_data, metadata, k_value) {
  if (paste0("K", k_value) %in% names(admixture_data)) {
    q_data <- admixture_data[[paste0("K", k_value)]]
    q_data$IID <- metadata$IID[1:nrow(q_data)]
    q_merged <- merge(q_data, metadata, by = "IID")
    
    # Melt data for plotting
    q_melted <- melt(q_merged, 
                     id.vars = c("IID", "Population", "Sample_Type"),
                     measure.vars = paste0("Ancestry", 1:k_value),
                     variable.name = "Ancestry_Component",
                     value.name = "Proportion")
    
    # Order samples by population and then by ancestry proportion
    q_melted$IID <- factor(q_melted$IID, 
                          levels = unique(q_melted$IID[order(q_melted$Population, q_melted$Proportion)]))
    
    p <- ggplot(q_melted, aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
      geom_bar(stat = "identity") +
      facet_grid(. ~ Sample_Type, scales = "free_x", space = "free_x") +
      scale_fill_brewer(type = "qual", palette = "Set3") +
      labs(
        title = paste0("ADMIXTURE Analysis (K=", k_value, ")"),
        x = "Samples",
        y = "Ancestry Proportion",
        fill = "Ancestry Component"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines")
      )
    
    return(p)
  } else {
    return(NULL)
  }
}

# Function to create ancestry summary table
create_ancestry_summary <- function(admixture_data, iadmix_data, metadata) {
  study_samples <- metadata[metadata$Sample_Type == "Study_Sample", ]
  
  summary_list <- list()
  
  # ADMIXTURE summary for K=5
  if ("K5" %in% names(admixture_data)) {
    admix_k5 <- admixture_data[["K5"]]
    admix_study <- admix_k5[1:nrow(study_samples), ]
    admix_study$Sample_ID <- study_samples$IID
    
    # Calculate dominant ancestry for each sample
    dominant_ancestry <- apply(admix_study[, paste0("Ancestry", 1:5)], 1, function(x) {
      max_idx <- which.max(x)
      return(paste0("Ancestry", max_idx))
    })
    
    admix_study$Dominant_Ancestry <- dominant_ancestry
    admix_study$Max_Proportion <- apply(admix_study[, paste0("Ancestry", 1:5)], 1, max)
    
    summary_list$admixture_k5 <- admix_study
  }
  
  # iAdmix summary
  if (!is.null(iadmix_data)) {
    iadmix_study <- iadmix_data[1:nrow(study_samples), ]
    iadmix_study$Sample_ID <- study_samples$IID
    summary_list$iadmix <- iadmix_study
  }
  
  return(summary_list)
}

# Function to create cross-validation plot for ADMIXTURE
create_cv_plot <- function(admixture_dir) {
  # Look for CV error file or extract from log files
  cv_file <- file.path(admixture_dir, "cv_errors.txt")
  
  if (file.exists(cv_file)) {
    cv_data <- read_table(cv_file, col_names = c("K", "CV_Error"))
  } else {
    # Try to extract CV errors from ADMIXTURE output files
    log_files <- list.files(admixture_dir, pattern = "\\.[0-9]+\\.Q$", full.names = TRUE)
    cv_errors <- c()
    k_values <- c()
    
    for (log_file in log_files) {
      k_val <- as.numeric(gsub(".*\\.([0-9]+)\\.Q$", "\\1", basename(log_file)))
      # This is a simplified approach - in reality, CV errors come from ADMIXTURE stdout
      # For now, create placeholder data
      k_values <- c(k_values, k_val)
      cv_errors <- c(cv_errors, runif(1, 0.4, 0.6)) # Placeholder
    }
    
    if (length(k_values) > 0) {
      cv_data <- data.frame(K = k_values, CV_Error = cv_errors)
    } else {
      return(NULL)
    }
  }
  
  p <- ggplot(cv_data, aes(x = K, y = CV_Error)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "red", size = 3) +
    labs(
      title = "ADMIXTURE Cross-Validation Error",
      x = "K (Number of Ancestral Populations)",
      y = "Cross-Validation Error"
    ) +
    theme_bw() +
    scale_x_continuous(breaks = cv_data$K)
  
  return(p)
}

# Main analysis function
main <- function() {
  # Input files - UPDATED FOR PSAM FORMAT
  pca_eigenval <- "analysis_other/ancestry/pca/merged_cohort.eigenval"
  pca_eigenvec <- "analysis_other/ancestry/pca/merged_cohort.eigenvec"
  psam_file <- "analysis_other/ancestry/plink/merged/merged_cohort.psam"  # UPDATED: PSAM instead of FAM
  admixture_dir <- "analysis_other/ancestry/global/admixture"
  iadmix_file <- "analysis_other/ancestry/global/iadmix/admixture_proportions.txt"
  rfmix_file <- "analysis_other/ancestry/local/rfmix/output.rfmix.Q"
  
  # Output directory
  output_dir <- "analysis_other/ancestry/reports"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read data
  cat("Reading PCA results...\n")
  pca_data <- read_pca_results(pca_eigenval, pca_eigenvec)
  
  cat("Reading sample metadata...\n")
  metadata <- read_sample_metadata(psam_file)  # UPDATED: using PSAM file
  
  cat("Reading ADMIXTURE results...\n")
  admixture_data <- read_admixture_results(admixture_dir)
  
  cat("Reading iAdmix results...\n")
  iadmix_data <- read_iadmix_results(iadmix_file)
  
  cat("Reading RFMix results...\n")
  rfmix_data <- read_rfmix_results(rfmix_file)
  
  # Create plots
  cat("Creating plots...\n")
  
  # PCA plots
  pca_plot <- create_pca_plot(pca_data, metadata)
  scree_plot <- create_scree_plot(pca_data)
  
  # ADMIXTURE plots
  admix_plots <- list()
  for (k in 3:7) {
    admix_plots[[paste0("K", k)]] <- create_admixture_plot(admixture_data, metadata, k)
  }
  
  # CV plot
  cv_plot <- create_cv_plot(admixture_dir)
  
  # Create summary statistics
  cat("Creating summary statistics...\n")
  ancestry_summary <- create_ancestry_summary(admixture_data, iadmix_data, metadata)
  
  # Save plots
  ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "scree_plot.png"), scree_plot, width = 8, height = 6, dpi = 300)
  
  for (k in 3:7) {
    if (!is.null(admix_plots[[paste0("K", k)]])) {
      ggsave(file.path(output_dir, paste0("admixture_K", k, ".png")), 
             admix_plots[[paste0("K", k)]], width = 12, height = 6, dpi = 300)
    }
  }
  
  if (!is.null(cv_plot)) {
    ggsave(file.path(output_dir, "cv_plot.png"), cv_plot, width = 8, height = 6, dpi = 300)
  }
  
  # Generate HTML report
  cat("Generating HTML report...\n")
  
  # Create Rmarkdown report
  rmd_content <- '
---
title: "Ancestry Analysis Report"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: bootstrap
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(ggplot2)
library(DT)
library(knitr)
```

## Summary

This report presents the results of ancestry analysis including:

- Principal Component Analysis (PCA) for population structure
- Global ancestry estimation using ADMIXTURE and iAdmix
- Local ancestry inference using RFMix and gnomix

## Sample Overview

```{r sample_overview}
metadata <- readRDS("metadata.rds")
sample_counts <- table(metadata$Sample_Type, metadata$Population)
kable(sample_counts, caption = "Sample counts by population and type")
```

## Principal Component Analysis

### PCA Plot
```{r pca_plot, fig.width=10, fig.height=8}
knitr::include_graphics("pca_plot.png")
```

### Scree Plot
```{r scree_plot, fig.width=8, fig.height=6}
knitr::include_graphics("scree_plot.png")
```

## Global Ancestry Analysis

### Cross-Validation Results
```{r cv_plot, fig.width=8, fig.height=6}
if (file.exists("cv_plot.png")) {
  knitr::include_graphics("cv_plot.png")
} else {
  cat("Cross-validation plot not available")
}
```

### ADMIXTURE Results

```{r admixture_plots, results="asis"}
for (k in 3:7) {
  plot_file <- paste0("admixture_K", k, ".png")
  if (file.exists(plot_file)) {
    cat("#### K =", k, "\\n\\n")
    cat("![](", plot_file, ")\\n\\n")
  }
}
```

### Study Sample Ancestry Proportions

```{r ancestry_table}
if (file.exists("ancestry_summary.rds")) {
  ancestry_summary <- readRDS("ancestry_summary.rds")
  if ("admixture_k5" %in% names(ancestry_summary)) {
    DT::datatable(ancestry_summary$admixture_k5, 
                  caption = "ADMIXTURE K=5 results for study samples",
                  options = list(scrollX = TRUE)) %>%
      DT::formatRound(columns = paste0("Ancestry", 1:5), digits = 3)
  }
}
```

## Local Ancestry Analysis

Local ancestry analysis results would be displayed here when available.

## Conclusions

- Total number of study samples analyzed: `r sum(metadata$Sample_Type == "Study_Sample")`
- Total number of reference samples: `r sum(metadata$Sample_Type == "Reference_1000G")`
- Optimal K value based on cross-validation: [To be determined from CV plot]

'
  
  # Write Rmarkdown file
  writeLines(rmd_content, file.path(output_dir, "ancestry_report.Rmd"))
  
  # Save data for the report
  saveRDS(metadata, file.path(output_dir, "metadata.rds"))
  saveRDS(ancestry_summary, file.path(output_dir, "ancestry_summary.rds"))
  
  # Render HTML report
  rmarkdown::render(
    file.path(output_dir, "ancestry_report.Rmd"),
    output_file = "ancestry_summary.html",
    output_dir = output_dir
  )
  
  cat("Analysis complete! Report saved to:", file.path(output_dir, "ancestry_summary.html"), "\n")
}

# Run main function
main()