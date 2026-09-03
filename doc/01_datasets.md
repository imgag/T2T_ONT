# Datasets used in this study

## Published sequencing data

We used publicly available sequencing data for benchmarking. Although assemblies are generated de novo without a reference genome, reference genomes are used to assess assembly quality.

Published datasets from Koren et al. 2024 (HG002):

- UL reads: https://s3-us-west-2.amazonaws.com/humanpangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/nanopore/ultra-long/
- Simplex reads: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/simplex/
- Pore-C: https://www.ncbi.nlm.nih.gov/sra/?term=SRR27664048

## Reference genomes

### Haploid reference genome (CHM13-T2T)

The latest release (T2T-CHM13v2.0) includes a Y chromosome and is used for mapping and comparison of final assemblies: [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)

### Diploid reference genome (HG002)

Published reference genome: [HG002 on GitHub](https://github.com/marbl/HG002?tab=readme-ov-file)

Current version: hg002v1.1 (July 2024)

Download: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz

Methods described in: McCartney et al. Nature Methods (June 2022): 687–95. https://doi.org/10.1038/s41592-022-01440-3

GIAB benchmark variant set (HG002 Q100 against CHM13v2.0):
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz

### Genome annotation

cDNA annotation for gene completeness assessment (Ensembl release 113):
https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/

Centromere / satellite annotation (CenSat v2.0):
https://github.com/hloucks/CenSatData

Downloaded from:
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.bed

Used for GQC error classification: `satellite` category (alpha-satellite HOR arrays, HSat1–3,
bSat, gSat, cenSat, ct, mixedAlpha, mon) and `rDNA` category (rDNA prefix entries covering
the 45S rDNA arrays on acrocentric chromosome short arms).

Gene annotation — JHU Liftoff v0.6 (CDS / rRNA on HG002v1.1):
Hansen et al., A complete diploid human genome benchmark for personalized genomics, bioRxiv 2025.
Two-pass Liftoff of RefSeq v110 from CHM13 onto HG002v1.1 (MAT and PAT separately).

Downloaded from:
```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/T2T/HG002/assemblies/annotation/JHULiftoff/v0.6/hg002.v1.1.loff.v0.6.mat.gff.gz \
  data/ref/hg002.v1.1.loff.v0.6.mat.gff.gz
aws s3 cp --no-sign-request \
  s3://human-pangenomics/T2T/HG002/assemblies/annotation/JHULiftoff/v0.6/hg002.v1.1.loff.v0.6.pat.gff.gz \
  data/ref/hg002.v1.1.loff.v0.6.pat.gff.gz
```

CDS features extracted from both haplotypes and merged into `data/ref/hg002v1.1.cds_jhu_liftoff.bed`
by the Snakemake rule `hg002_cds_bed`. Used for GQC error classification `coding` category.



Used for GQC error classification `segmental_dup` category.

Tandem repeat annotation (GQC resource files — MNR ≥10 bp, DNR, TNR, QNR):
Pre-computed by the GQC pipeline as part of `data/ref/gqc_resources/`.
Files: `v1.1.mnr10.bed`, `v1.1.dnr.bed`, `v1.1.tnr.bed`, `v1.1.qnr.bed`.
Merged by the Snakemake rule `hg002_tandem_repeat_bed` into `data/ref/hg002v1.1.tandem_repeats.bed`.
Used for GQC error classification `tandem_repeat` category.

## Repeat Masker database

DFAM database for mammals/human in H5 format, downloaded from https://www.dfam.org/releases/current/families/FamDB/

## Segmental Duplications track

### CHM13:

From Vollger et al. 2022, Cell:
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.full.bed

### HG002
Segmental duplications (SEDEF, hg002v1.1, Jan 2025):

Downloaded from:
```bash
aws s3 cp --no-sign-request \
  s3://human-pangenomics/T2T/HG002/assemblies/annotation/segdups/hg002v1.1.SDs.013025.bed \
  data/ref/hg002v1.1.SDs.013025.bed
```

## HPRC comparison assemblies

The metadata file was downloaded from `data_tables/assemblies_pre_release_v0.6.1.index.csv` from Feb 14, 2025 (Commit #41aa47d)  from the https://github.com/human-pangenomics/hprc_intermediate_assembly GitHub.

https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/41aa47dd3430fbb250cdb6a78efde43313d35557/data_tables/assemblies_pre_release_v0.6.1.index.csv


HPRC Release 2 assemblies were downloaded for comparison using the provided script:

```bash
python workflow/scripts/22_download_hprc_asm.py \
    --output data/ref/hprc \
    --qc-output-dir assembly/qc/phased_hprc \
    --config-file data/ref/hprc/download_status.txt \
    --flagger-file data/ref/hprc/flagger_ont_v0.1.1.csv \
    --nucflag-file data/ref/hprc/nucflag_pre_release_v0.1.index.csv \
    doc/tables/HPGRC_R2_SequencingDataIndex_assemblies.tsv 
```

The following 20 HPRC samples were randomly selected and downloaded (note: HG002 was excluded from assembly quality comparison plots, as it is also the diploid reference genome and):

```
HG002   HG00272  HG005    HG00639  HG01252  HG01361  HG02083  HG02155
HG03834 HG04228  HG06807  NA18522  NA18747  NA18948  NA19043  NA20503
NA20762 NA20806  NA20809  NA20827
```

## HG002 low complexity regions

Bedfiles for TandemRepeats and Homopolymers published by NIST. We used v3.6 with separate files for the paternal and maternal haplotype of HG002.




## 1000 Genomes variant sets

Phased 1000 Genomes variants (Lalli et al.), mapped to CHM13v2.0:
https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1

## dbSNP variant sets

### CHM13-T2T

CHM13v2.0 dbSNP build 155, lifted over from GRCh38:
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_dbSNPv155.vcf.gz
