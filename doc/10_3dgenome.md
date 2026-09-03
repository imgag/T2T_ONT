# 3D Genome Analysis

Documentation of haplotype-resolved 3D genome reconstruction from Pore-C data.

## Haplotype-specific chromatin contacts

Pore-C reads are haplotagged using phased variant calls (WhatsHap), enabling separation of contacts into maternal and paternal haplotypes. This allows detection of allele-specific chromatin organization, including X-chromosome inactivation (XCI)-related structures.

Notable allele-specific features identified in the literature (Chen et al. 2025):

- X chromosome inactivation superdomain at ~115 Mb, with a superloop between tandem repeats DXZ4 and FIRRE
- Allele-specific sub-TAD organization at imprinted loci (Dlk1-Dio3, Igf2-H19)

## Metrics

- **Stratum-adjusted correlation coefficient (SCC)**: correlation between Pore-C and Hi-C contact matrices, computed using HiCrep
- **Compartmental eigenvector scores**: computed with `cooltools eigs-cis`
- **TAD insulation scores**: computed with `cooltools diamond-insulation`

## Implementation

3D genome analysis is implemented in:

- `workflow/rules/dip3d.smk`: Dip3D haplotype-resolved contact pipeline
- `workflow/rules/porec.smk`: Pore-C alignment and pairs generation
- `workflow/rules/structure3d.smk`: 3D structure prediction (ASHIC, LorDG)
