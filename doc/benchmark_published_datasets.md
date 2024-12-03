# Benchmark on published datasets

## Raw data:

From Koren et al 2024:
- UL:  https://s3-us-west-2.amazonaws.com/humanpangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/nanopore/ultra-long/
- Duplex: https://s3-us-west2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/
- Pore-C: https://www.ncbi.nlm.nih.gov/sra/?term=SRR27664048
- Simplex: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/simplex/

## Comparisons

Use benchmark data to compare approaches with different HQ datasets:

### HQ Datasets

- Comparison HERRO to Duplex

### UL datasets

- Our flowcells to published dataset
- Yield, Qscore, N50


## Pipeline validations

### Subsampling comparison UL

- Does coverage fit the expected?
- How does N50 and QScore distribution improve with subsampling?


### ONT T2T-Sequencing Kit

ONT has published a protocol for T2T assembly on their protocols page [here](https://nanoporetech.com/document/telomere-to-telomere-sequencing). However this does not include detailed analyses steps, they mention a complete workflow will be published soon.

