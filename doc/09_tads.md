# TAD and Loop Analysis

Documentation of 3D chromatin structure analysis using Pore-C data.

## Overview

Haplotype-resolved chromatin contacts are generated from Pore-C sequencing. Each haplotype is processed separately to produce phased contact matrices. Analysis is implemented in `workflow/rules/dip3d.smk` and `workflow/rules/porec.smk`.

## Definitions

**TAD (Topologically Associating Domain) calling:**

TADs are large genomic regions (typically 100 kb–1 Mb) where DNA sequences interact more frequently with each other than with sequences outside the domain. They represent stable chromatin organization units conserved across cell types. TAD boundaries are identified using insulation score calculations.

**Loop calling:**

Loops are specific point-to-point contacts between distant genomic loci (e.g., enhancer–promoter interactions). They represent dynamic, often tissue-specific regulatory interactions and are identified as statistically significant peak interactions above local background.

## Analysis approach

1. Pore-C reads are aligned and processed using `wf-pore-c` or `pairtools`.
2. Contact matrices are generated in `.pairs` and `.cool` formats.
3. Matrices are balanced using `cooler balance`.
4. Compartment analysis: `cooltools eigs-cis` (100 kb resolution).
5. Insulation score and TAD boundary calling: `cooltools diamond-insulation`.
6. Quality metrics: stratum-adjusted correlation coefficient using [HiCrep](https://github.com/TaoYang-dev/hicrep).

## TAD callers

Tools evaluated for TAD calling:
- [cooltools](https://github.com/open2c/cooltools) (primary)
- OnTAD
- Armatus
- Arrowhead (Juicer Tools)

## Additional tools

- [dCHiC](https://github.com/ay-lab/dcHiC): differential Hi-C analysis
- [HiCrep](https://github.com/TaoYang-dev/hicrep): stratum-adjusted correlation between contact matrices
