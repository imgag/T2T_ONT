# R presets for plotting
# Source this file at the beginning of your R scripts

library(ggplot2)
suppressPackageStartupMessages(library(extrafont))

# Recommendation (Nature Guidelines)
# Font: Helvetica (substituted here with Liberation Sans, a free
# metric-compatible equivalent -- see doc/figures/helvetica/README.md)
# Multi Part Figure Labels: 8pt bold, a b c d ..
# Max text size: 7pt
# Min text size: 5pt

font_import(paths = "../../doc/figures/helvetica", prompt = FALSE)
loadfonts(device = "pdf") # or device = "win" for Windows

# Custom ggplot theme based on Nature guidelines
# Usage: theme_set(t2t_theme()) or customize with t2t_theme(base_size = 7)
t2t_theme <- function(base_size = 6, base_family = "Liberation Sans") {
    theme_bw(base_size = base_size, base_family = base_family) +
        theme(
            # Text elements (in pt)
            plot.title = element_text(size = 8, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 7, hjust = 0),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6),
            legend.title = element_text(size = 7, face = "bold"),
            legend.text = element_text(size = 6),
            strip.text = element_text(size = 7, face = "bold"),
            
            # Background and grid
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.background = element_blank(),
            panel.spacing = unit(0.5, "lines"),
            
            # Plot margins - add padding around data points
            plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt"),
            
            # Legend
            legend.position = "right",
            legend.key.size = unit(0.8, "lines")
        )
}

# Set as default theme
theme_set(t2t_theme())

# Define a color palettes
col_source <- c("HPRC" = "#9ba39f", "Samples" = "#bb1133")


col_hp1 <- "#4292c6"
col_hp2 <- "#ef6548"
col_diploid <- "#787878"

col_haplotype <- c(
    "Hap 1" = col_hp1,
    "Hap 2" = col_hp2,
    "haplotype1" = col_hp1,
    "haplotype2" = col_hp2,
    "paternal" = col_hp1,
    "maternal" = col_hp2,
    "hap1" = col_hp1,
    "hap2" = col_hp2,
    "both" = col_diploid,
    "diploid" = col_diploid,
    "h1_frags" = col_hp1,
    "h2_frags" = col_hp2,
    "hp1" = col_hp1,
    "hp2" = col_hp2,
    "unphased_frags" = col_diploid
)


superpop_colors <- c(
  "EUR" = "#1976D2",  # Blue
  "AFR" = "#FF6F00",  # Orange
  "EAS" = "#388E3C",  # Green
  "SAS" = "#D32F2F",  # Red
  "AMR" = "#7B1FA2"   # Purple
)

assign_superpop <- function(country) {
  if (country %in% c("Germany", "Austria", "France", "Spain", "Croatia", "Poland", "Russia")) {
    return("EUR")
  } else if (country %in% c("Taiwan", "Vietnam")) {
    return("EAS")
  } else if (country %in% c("India")) {
    return("SAS")
  } else if (country %in% c("Mexico", "Chile")) {
    return("AMR")
  } else {
    return("EUR")
  }
}



# Create karyotype data from reference FAI file
chm13 <- "../../data/ref/T2T-CHM13.v2.fasta.fai"
karyotype_chm13 <- read_tsv(chm13, 
    col_names = c("chr", "end", "offset", "linebases", "linewidth", "qualoffset"),
    col_types = "cnnnnn",
    show_col_types = FALSE) %>%
    filter(chr != "chrM") %>%
    mutate(start = 1) %>%
    select(chr, start, end) %>%
    # Create numeric sort key BEFORE any factor conversion
    mutate(chr_num = case_when(
        chr == "chrX" ~ 23,
        chr == "chrY" ~ 24,
        TRUE ~ as.numeric(str_remove(chr, "chr"))
    )) %>%
     mutate(
        chr = factor(chr, levels = unique(chr[order(chr_num)])),
    ) %>%
    arrange(chr_num)

# Duplicate haplotypes for diploid plotting
karyotype_chm13_diploid <- karyotype_chm13 %>%
    crossing(haplotype = c("Hap 1", "Hap 2")) %>%
    # Now create ordered factors
    mutate(
        chr = factor(chr, levels = unique(chr[order(chr_num)])),
        haplotype = factor(haplotype, levels = c("Hap 1", "Hap 2"))
    ) %>%
    select(-chr_num) %>%
    arrange(chr, haplotype) %>%
    filter(!(chr=="chrY" & haplotype=="Hap 2")) # Keep chrY only for Hap 1

# Read cytoband file
cytoband <- read_tsv("../../data/ref/cytoBandMapped.bed",
                     col_names = c("chr", "start", "end", "name", "gieStain", "full_name"),
                     show_col_types = FALSE) %>%
            separate_wider_position(name, 
                widths = c("arm" = 1),
                too_many = "drop",
                cols_remove = FALSE ) 


qc_table_to_wide <- function(dt, col_order = NULL) {
    dt_qc_table <- dt %>%
    filter() %>%
    mutate(haplotype = ifelse(haplotype == "assembly", "both", haplotype)) %>%
    group_by(asm_name, metric, haplotype) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(
        names_from = metric, 
        values_from = value  # Keep as numeric, don't format here
    )
    
    # Reorder columns if col_order is provided
    if (!is.null(col_order)) {
        # Keep asm_name and haplotype first, then reorder metric columns
        existing_cols <- intersect(col_order, colnames(dt_qc_table))
        dt_qc_table <- dt_qc_table %>%
            select(asm_name, haplotype, all_of(existing_cols))
    }
    
    return(dt_qc_table)
}

metrics_phased <- c(
    "length",
    "n_contigs_over_10mb",
    "n_contigs",
    "ref_covered",
    "asm_covered",
    "genome_completeness",
    "MMC",
    "MSC",
    "intersection_blocks",
    "all_switch_rate",
    "blockwise_hamming_rate",
    "n_t2t_contig",
    "n_t2t_scaffold",
    "total_n_count",
    "total_gaps"
)

metrics_full <- c(
    "length",
    "n_contigs_over_10mb",
    "n_contigs",
    "n_t2t_contig",
    "n_t2t_scaffold",
    "ref_covered",
    "asm_covered",
    "genome_completeness",
    "MMC",
    "MSC",
    "inter-chromosomal_misjoins",
    "intra-chromosomal_gaps",
    "candidate_inversions_in_the_middle",
    "candidate_inversions_at_contig_ends",
    "total_gaps",
    "total_n_count",
    "total_gaps",
    "total_covered_variants",
    "overall_switch_rate",
    "overall_switch_flip_rate",
    "total_blocks",
    "total_covered_variants",
    "total_switches",
    "total_het_variants",
    "overall_switch_rate",
    "overall_switchflip_rate",
    "overall_hamming_rate"
)

# Metrics summary is condensed
# Only values for two haplotypes, not both.
metrics_summary <- c(
    "length",
    "n_t2t_contig",
    "n_t2t_scaffold",
    "n_contigs_over_10mb",
    "genome_completeness",
    "MSC",
    "MMC",
    "ref_covered",
    "asm_covered"
)


# Shared column name lookup for gt (HTML) and Excel export
# Used by both qc_table_pretty() and export_excel()
col_labels <- list(
  # GT format (with HTML line breaks)
  gt = c(
    "asm_name" = "Assembly",
    "haplotype" = "Haplotype",
    "length" = "Assembly<br>Length",
    "n_contigs" = "n<br>Contigs",
    "n_contigs_over_10mb" = "Contigs<br>>10Mbp",
    "genome_completeness" = "Genome<br>Completeness",
    "MSC" = "MSC",
    "MMC" = "MMC",
    "ref_covered" = "Reference<br>Covered",
    "asm_covered" = "Assembly<br>Covered",
    "n_t2t_contig" = "T2T<br>Contigs",
    "n_t2t_scaffold" = "T2T<br>Scaffolds",
    "inter-chromosomal_misjoins" = "Trans-chr<br>Misjoins",
    "intra-chromosomal_gaps" = "Cis-chr<br>Gaps",
    "candidate_inversions_in_the_middle" = "Inversions<br>Middle",
    "candidate_inversions_at_contig_ends" = "Inversions<br>Ends",
    "total_gaps" = "Total<br>Gaps",
    "total_n_count" = "N<br>Count",
    "total_covered_variants" = "Covered<br>Variants",
    "overall_switch_rate" = "Switch<br>Rate",
    "total_blocks" = "Phase<br>Blocks",
    "total_switches" = "Switches",
    "total_het_variants" = "Het<br>Variants",
    "overall_switchflip_rate" = "Switchflip<br>Rate",
    "overall_hamming_rate" = "Hamming<br>Rate",
    "overall_diff_genotypes_rate" = "Diff Genotypes<br>Rate",
    "intersection_blocks" = "Intersection<br>Blocks",
    "all_switch_rate" = "All Switch<br>Rate",
    "blockwise_hamming_rate" = "Blockwise<br>Hamming Rate"
  ),
  # Excel format (with newline breaks)
  excel = c(
    "asm_name" = "Assembly\nName",
    "haplotype" = "Haplotype",
    "length" = "Assembly\nLength",
    "n_contigs" = "n\nContigs",
    "n_contigs_over_10mb" = "Contigs\n>10Mbp",
    "genome_completeness" = "Genome\nCompleteness",
    "MSC" = "MSC",
    "MMC" = "MMC",
    "ref_covered" = "Reference\nCovered",
    "asm_covered" = "Assembly\nCovered",
    "n_t2t_contig" = "T2T\nContigs",
    "n_t2t_scaffold" = "T2T\nScaffolds",
    "inter-chromosomal_misjoins" = "Trans-chr\nMisjoins",
    "intra-chromosomal_gaps" = "Cis-chr\nGaps",
    "candidate_inversions_in_the_middle" = "Inversions\nMiddle",
    "candidate_inversions_at_contig_ends" = "Inversions\nEnds",
    "total_gaps" = "Total\nGaps",
    "total_n_count" = "N\nCount",
    "total_covered_variants" = "Covered\nVariants",
    "overall_switch_rate" = "Switch\nRate",
    "total_blocks" = "Phase\nBlocks",
    "total_switches" = "Switches",
    "total_het_variants" = "Het\nVariants",
    "overall_switchflip_rate" = "Switchflip\nRate",
    "overall_hamming_rate" = "Hamming\nRate",
    "overall_diff_genotypes_rate" = "Diff\nGenotypes\nRate",
    "intersection_blocks" = "Intersection\nBlocks",
    "all_switch_rate" = "All Switch\nRate",
    "blockwise_hamming_rate" = "Blockwise\nHamming Rate"
  )
)

# Column formatting rules
col_format <- list(
  # Columns to format as Gbp (scale by 1e-9, 2 decimals)
  gbp = c("length"),
  # Columns to format as fractions (3 decimals)
  fraction = c("genome_completeness", "ref_covered", "asm_covered", 
               "overall_switch_rate", "overall_switchflip_rate", 
               "overall_hamming_rate", "overall_diff_genotypes_rate",
               "all_switch_rate", "blockwise_hamming_rate"),
  # Columns to format as decimals (2 decimals)
  decimal = c("MMC", "MSC"),
  # Columns to format as integers (0 decimals)

  integer = c("n_contigs", "n_contigs_over_10mb", "n_t2t_contig", "n_t2t_scaffold",
              "total_gaps", "total_n_count", "inter-chromosomal_misjoins",
              "intra-chromosomal_gaps", "candidate_inversions_in_the_middle",
              "candidate_inversions_at_contig_ends", "total_covered_variants",
              "total_blocks", "total_switches", "total_het_variants",
              "intersection_blocks")
)

qc_table_pretty <- function(qc_selected, title = "Assembly QC Metrics", 
                            subtitle = "Selected metrics") {
    
  # Get column names present in the data (excluding asm_name and haplotype)
  data_cols <- setdiff(colnames(qc_selected), c("asm_name", "haplotype"))
  
  # Build dynamic column labels for present columns
  gt_labels <- col_labels$gt[intersect(names(col_labels$gt), colnames(qc_selected))]
  
  # Prepare data
  qc_prepared <- qc_selected %>%
    arrange(asm_name, haplotype) %>%
    group_by(asm_name) %>%
    mutate(asm_name = ifelse(row_number() == 1, asm_name, "")) %>%
    ungroup()
  
  # Start building the gt table
  qc_gt_table <- qc_prepared %>%
    gt(rowname_col = "haplotype") %>%
    tab_header(
      title = md(paste0("**", title, "**")),
      subtitle = subtitle
    )
  
  # Apply column labels dynamically
  for (col_name in names(gt_labels)) {
    if (col_name %in% colnames(qc_selected)) {
      qc_gt_table <- qc_gt_table %>%
        cols_label(!!sym(col_name) := html(gt_labels[col_name]))
    }
  }
  
  # Apply formatting for Gbp columns
  gbp_cols <- intersect(col_format$gbp, data_cols)
  if (length(gbp_cols) > 0) {
    qc_gt_table <- qc_gt_table %>%
      fmt_number(
        columns = all_of(gbp_cols),
        decimals = 2,
        scale_by = 1/1e9,
        pattern = "{x} Gbp"
      )
  }
  
  # Apply formatting for fraction columns
  frac_cols <- intersect(col_format$fraction, data_cols)
  if (length(frac_cols) > 0) {
    qc_gt_table <- qc_gt_table %>%
      fmt_number(
        columns = all_of(frac_cols),
        decimals = 3
      )
  }
  
  # Apply formatting for decimal columns
  dec_cols <- intersect(col_format$decimal, data_cols)
  if (length(dec_cols) > 0) {
    qc_gt_table <- qc_gt_table %>%
      fmt_number(
        columns = all_of(dec_cols),
        decimals = 2
      )
  }
  
  # Apply formatting for integer columns
  int_cols <- intersect(col_format$integer, data_cols)
  if (length(int_cols) > 0) {
    qc_gt_table <- qc_gt_table %>%
      fmt_number(
        columns = all_of(int_cols),
        decimals = 0
      )
  }
  
  # Apply common styling
  qc_gt_table <- qc_gt_table %>%
    fmt_missing(
      columns = everything(),
      missing_text = "—"
    ) %>%
    cols_align(
      align = "center",
      columns = all_of(data_cols)
    ) %>%
    cols_align(
      align = "left",
      columns = asm_name
    ) %>%
    tab_options(
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      table.width = pct(90),
      heading.title.font.size = px(20),
      heading.subtitle.font.size = px(15),
      column_labels.border.bottom.color = "black",
      column_labels.font.weight = "bold",
      table_body.hlines.color = "#D3D3D3"
    ) %>%
    # Add thick top border to first row of each group
    tab_style(
      style = cell_borders(
        sides = "top",
        color = "black",
        weight = px(4)
      ),
      locations = cells_body(
        rows = qc_selected %>%
          arrange(asm_name, haplotype) %>%
          group_by(asm_name) %>%
          mutate(is_first = row_number() == 1) %>%
          ungroup() %>%
          pull(is_first)
      )
    )
  
  return(qc_gt_table)
}


export_excel <- function(dt, outfile) {

  # Export full table to Excel with formatting
  library(openxlsx)
  
  # Get column names present in the data
  data_cols <- colnames(dt)

  # Prepare formatted data for Excel
  qc_excel <- dt
  
  # Apply Gbp formatting
  gbp_cols <- intersect(col_format$gbp, data_cols)
  for (col in gbp_cols) {
    qc_excel[[col]] <- round(as.numeric(qc_excel[[col]]) / 1e9, 2)
  }
  
  # Apply fraction formatting
  frac_cols <- intersect(col_format$fraction, data_cols)
  for (col in frac_cols) {
    qc_excel[[col]] <- round(as.numeric(qc_excel[[col]]), 3)
  }
  
  # Apply decimal formatting
  dec_cols <- intersect(col_format$decimal, data_cols)
  for (col in dec_cols) {
    qc_excel[[col]] <- round(as.numeric(qc_excel[[col]]), 2)
  }
  
  # Apply integer formatting
  int_cols <- intersect(col_format$integer, data_cols)
  for (col in int_cols) {
    qc_excel[[col]] <- as.integer(qc_excel[[col]])
  }

  # Create workbook
  wb <- createWorkbook()
  addWorksheet(wb, "QC Table")

  # Write data
  writeData(wb, "QC Table", qc_excel)

  # Style header row with text wrapping
  headerStyle <- createStyle(
    fontName = "Helvetica",
    fontSize = 11,
    fontColour = "black",
    halign = "center",
    valign = "center",
    textDecoration = "bold",
    border = "bottom",
    borderColour = "black",
    wrapText = TRUE
  )
  addStyle(wb, "QC Table", headerStyle, rows = 1, cols = 1:ncol(qc_excel), gridExpand = TRUE)

  # Set header row height to accommodate wrapped text
  setRowHeights(wb, "QC Table", rows = 1, heights = 40)

  # Update header names using shared lookup
  for (i in seq_along(colnames(qc_excel))) {
    col_name <- colnames(qc_excel)[i]
    if (col_name %in% names(col_labels$excel)) {
      writeData(wb, "QC Table", col_labels$excel[col_name], 
                startCol = i, startRow = 1)
    }
  }

  # Style data columns
  dataStyle <- createStyle(
    fontName = "Helvetica",
    fontSize = 10,
    halign = "center"
  )
  addStyle(wb, "QC Table", dataStyle, rows = 2:(nrow(qc_excel) + 1), 
          cols = 1:ncol(qc_excel), gridExpand = TRUE)

  # Format specific columns with units (Gbp)
  for (col in gbp_cols) {
    col_idx <- which(colnames(qc_excel) == col)
    if (length(col_idx) > 0) {
      addStyle(wb, "QC Table", 
              createStyle(numFmt = "0.00 \"Gbp\""),
              rows = 2:(nrow(qc_excel) + 1), cols = col_idx, stack = TRUE)
    }
  }

  # Format fraction columns (3 decimals)
  frac_col_idx <- which(colnames(qc_excel) %in% frac_cols)
  if (length(frac_col_idx) > 0) {
    addStyle(wb, "QC Table",
            createStyle(numFmt = "0.000"),
            rows = 2:(nrow(qc_excel) + 1), cols = frac_col_idx, gridExpand = TRUE, stack = TRUE)
  }

  # Format decimal columns (MMC, MSC)
  dec_col_idx <- which(colnames(qc_excel) %in% dec_cols)
  if (length(dec_col_idx) > 0) {
    addStyle(wb, "QC Table",
            createStyle(numFmt = "0.00"),
            rows = 2:(nrow(qc_excel) + 1), cols = dec_col_idx, gridExpand = TRUE, stack = TRUE)
  }

  # Format integer columns
  int_col_idx <- which(colnames(qc_excel) %in% int_cols)
  if (length(int_col_idx) > 0) {
    addStyle(wb, "QC Table",
            createStyle(numFmt = "0"),
            rows = 2:(nrow(qc_excel) + 1), cols = int_col_idx, gridExpand = TRUE, stack = TRUE)
  }

  # Add borders around groups of 3 rows (same asm_name)
  n_assemblies <- nrow(qc_excel) / 3
  if (n_assemblies >= 1 && n_assemblies == floor(n_assemblies)) {
    for (i in 1:n_assemblies) {
      start_row <- 2 + (i - 1) * 3
      end_row <- start_row + 2
      
      # Top border
      addStyle(wb, "QC Table",
              createStyle(border = "top", borderColour = "black", borderStyle = "medium"),
              rows = start_row, cols = 1:ncol(qc_excel), gridExpand = TRUE, stack = TRUE)
      
      # Bottom border
      addStyle(wb, "QC Table",
              createStyle(border = "bottom", borderColour = "black", borderStyle = "medium"),
              rows = end_row, cols = 1:ncol(qc_excel), gridExpand = TRUE, stack = TRUE)
      
      # Left border
      addStyle(wb, "QC Table",
              createStyle(border = "left", borderColour = "black", borderStyle = "medium"),
              rows = start_row:end_row, cols = 1, gridExpand = TRUE, stack = TRUE)
      
      # Right border
      addStyle(wb, "QC Table",
              createStyle(border = "right", borderColour = "black", borderStyle = "medium"),
              rows = start_row:end_row, cols = ncol(qc_excel), gridExpand = TRUE, stack = TRUE)
    }
  }

  # Set column widths
  setColWidths(wb, "QC Table", cols = 1:ncol(qc_excel), widths = "auto")

  # Save workbook
  saveWorkbook(wb, outfile, overwrite = TRUE)
}