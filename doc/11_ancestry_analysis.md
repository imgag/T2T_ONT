# Ancestry Analysis

Documentation of population genetics and ancestry inference methods used in this study.

## SNP requirements by method

### Global ancestry analysis

| **Method**    | **Typical SNPs** | **Range** | **Examples from literature**                              |
| ------------- | ---------------- | --------- | --------------------------------------------------------- |
| **ADMIXTURE** | 100K–1M          | 50K–2M    | 1000 Genomes: ~650K SNPs; UK Biobank: ~700K SNPs         |
| **STRUCTURE** | 50K–500K         | 10K–1M    | HapMap studies: ~300K SNPs                               |
| **iAdmix**    | 50K–300K         | 20K–500K  | Original paper: ~200K SNPs                               |
| **PCA-based** | 100K–1M          | 50K–2M    | EIGENSTRAT: ~500K SNPs                                   |

### Local ancestry analysis

| **Method** | **Typical SNPs** | **Range** | **Notes**                     |
| ---------- | ---------------- | --------- | ----------------------------- |
| **RFMix**  | 500K–5M          | 100K–10M  | Requires dense coverage       |
| **LAMP**   | 300K–2M          | 100K–5M   | Sliding window approach       |
| **Gnomix** | 1M–10M           | 500K–20M  | Deep learning approach        |

## Reference datasets

Landmark studies used for comparison:

1. **1000 Genomes Project (2015)**: ~650,000 high-quality SNPs for global ancestry inference
2. **Bryc et al. (2015)**: ~750,000 SNPs, ADMIXTURE-based global ancestry in African Americans
3. **Moreno-Estrada et al. (2013)**: ~300,000 SNPs, ADMIXTURE and RFMix for Latin American populations
4. **Hellenthal et al. (2014)**: ~650,000 SNPs, ChromoPainter/GLOBETROTTER across 95 populations worldwide
5. **UK Biobank (2018–present)**: ~700,000 genotyped SNPs

## Implementation

Ancestry analysis is implemented in `workflow/rules/ancestry.smk`. The pipeline uses phased VCF files aligned to CHM13v2.0 and the phased 1000 Genomes reference panel (Lalli et al., CHM13v2.0 coordinates).
