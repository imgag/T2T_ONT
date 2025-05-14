# RepeatMasker Analysis

This directory contains scripts to run and analyze RepeatMasker output.

## Contents

- `repeatmasker.sh`: Script to run RepeatMasker on an assembly
- `TUE_02_03UL/`: Directory containing RepeatMasker output for the TUE_02_03UL assembly

## RepeatMasker Output Files

RepeatMasker generates several output files:
- `.out`: Main output file with detailed annotation of each repeat
- `.tbl`: Summary table of repeat content
- `.out.gff`: Repeat annotations in GFF format
- `.masked`: Masked version of the input sequence
- `.out.html`: HTML report
- `.cat.gz`: Categorized output

## Using the Analysis Script

The analysis script is located at `workflow/scripts/16_analyze_repeatmasker.R`. It processes RepeatMasker `.out` files to:
1. Summarize repeat content by class and family
2. Generate visualizations of repeat distributions
3. Create BED files for downstream analysis (both filtered and unfiltered)

### Prerequisites

Required R packages:
```
tidyverse
ggplot2
viridis
patchwork
optparse
```

Install them with:
```R
install.packages(c("tidyverse", "ggplot2", "viridis", "patchwork", "optparse"))
```

### Usage

```bash
Rscript workflow/scripts/16_analyze_repeatmasker.R --input <path_to_out_file> [options]
```

### Options

- `-i, --input`: Input RepeatMasker .out file (required)
- `-o, --output`: Output prefix for generated files (default: "repeatmasker_analysis")
- `-m, --min_length`: Minimum length of repeats to include (default: 0)
- `-d, --max_div`: Maximum divergence percentage to include (default: 100)

### Example

```bash
Rscript workflow/scripts/16_analyze_repeatmasker.R --input TUE_02_03UL/assembly.fasta.out --output TUE_02_03UL_analysis --min_length 100 --max_div 25
```

### Output Files

The script generates the following output files:

#### Summary Files
- `<prefix>_class_summary.csv`: Summary statistics by repeat class
- `<prefix>_family_summary.csv`: Summary statistics by repeat family
- `<prefix>_sequence_summary.csv`: Summary statistics by sequence/chromosome

#### BED Files
- `<prefix>_unfiltered.bed`: BED file with all repeat annotations
- `<prefix>_unfiltered_simplified.bed`: Simplified BED file with all repeats
- `<prefix>_filtered.bed`: BED file with filtered repeat annotations
- `<prefix>_filtered_simplified.bed`: Simplified BED file with filtered repeats
- `<prefix>_nucplot.bed`: BED file formatted for nucplot visualization (filtered data)

#### Visualization Files
- `<prefix>_plots.pdf/png`: Combined visualization of main plots
- `<prefix>_class_summary.pdf/png`: Distribution of repeat classes
- `<prefix>_family_summary.pdf/png`: Top 20 repeat families by length
- `<prefix>_divergence_distribution.pdf/png`: Divergence distribution by repeat class
- `<prefix>_length_distribution.pdf/png`: Length distribution by repeat class (boxplot)
- `<prefix>_length_histogram.pdf/png`: Length histogram by repeat class
- `<prefix>_divergence_vs_length.pdf/png`: Divergence vs length scatter plot

## Visualizations

The script generates several plots to analyze repeat distributions:

### Basic Plots
1. Distribution of repeat classes by total length
2. Top 20 repeat families by total length
3. Divergence distribution by repeat class
4. Length distribution by repeat class (log scale boxplot)

### Additional Plots
5. Length histogram by repeat class (faceted)
6. Divergence vs length scatter plot (faceted)

## BED File Formats

### Standard BED Files
The standard BED files include the following columns:
1. Sequence name (chromosome)
2. Start position (0-based)
3. End position
4. Repeat name (class:family)
5. Repeat name (specific repeat name)
6. Strand
7. Start position (for thickness)
8. End position (for thickness)
9. RGB color (based on repeat class)

### Simplified BED Files
The simplified BED files include:
1. Sequence name (chromosome)
2. Start position (0-based)
3. End position
4. Repeat class

### Nucplot BED File
The nucplot BED file is specifically formatted for use with nucplot visualization tools:
1. Sequence name (chromosome)
2. Start position (0-based)
3. End position
4. Repeat class
5. Plot flag (always set to "plot")

Example:
```
sequence1  1000  1100  LINE  plot
sequence1  1200  1300  SINE  plot
```

## Notes

- Processing large RepeatMasker output files can require significant memory
- For very large files, consider filtering by sequence/chromosome first
- The script is designed to work with RepeatMasker output format version 4.1.0
- Color coding in BED files:
  - SINE: Red (255,0,0)
  - LINE: Green (0,255,0)
  - LTR: Blue (0,0,255)
  - DNA: Yellow (255,255,0)
  - Simple_repeat: Magenta (255,0,255)
  - Satellite: Cyan (0,255,255)
  - Low_complexity: Gray (128,128,128)