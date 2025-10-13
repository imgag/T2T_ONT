library(tidyverse)

f <- "/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/qc/phased_verkko/TUE_02_03UL/haplotype1.mapped_T2T.paf"

# PAF file column definitions
colnames <- c(
  "query_name",        # Query sequence name
  "query_length",      # Query sequence length
  "query_start",       # Query start (0-based; BED-like; closed)
  "query_end",         # Query end (0-based; BED-like; open)
  "strand",            # Relative strand: "+" or "-"
  "target_name",       # Target sequence name
  "target_length",     # Target sequence length
  "target_start",      # Target start on original strand (0-based)
  "target_end",        # Target end on original strand (0-based)
  "num_matches",       # Number of residue matches
  "aln_block_length",  # Alignment block length
  "mapping_quality"    # Mapping quality (0-255; 255 for missing)
)

paf <- read_tsv(f, col_names = colnames)

paf %>%
  filter(target_name == "chr13") %>%
  arrange(num_matches)

