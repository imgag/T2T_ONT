# Haplotype-phased T2T-3D genomes and genomewide topological domains of 20 individuals

```
Title comments:

- Maybe we should include near-T2T or similar on the title.
- We should include nanopore-only
- 3D Genome is good, topological domains are not relevant enough for title

Suggestion:

Near complete 3D-T2T genomes using Nanopore only for 20 individuals.
```

## Introduction

see external markdown doc

## Methods

### Cohort

Describe all participants and samples
In total we analysed the entire genomes of ??? adult individuals from different nationalities. Probands were from different European countries as Germany (???), Russia (?), Austria (?), Spain (?), and France (?), but also from countries outside Europe like Chile (?), Mexico (?), South India (?), and East Asia (???). Several individuals had known mixed ancestries (like Polish/Croatia/Serbia/Hungary; or German/Austrian, Eastslavic/Mordvine). Of the ??? individuals, ??? were females, the other ??? were males. Mean age was ?? with the range of ??? to … All probands were individually informed on the aims of the project and signed an agreement to openly share their genome data.  
From all  individuals ??ml EDTA and ??? heparinized blood were taken.

### DNA and PBMC preparation

### PoreC protocol

### Library prep and sequencing

Describe all Nanopore libraries, sequencing depth etc.

### Nanopore-only T2T genome assemblies

Library types, error correction, verkko assembly, poreC phasing

### Benchmarking T2T assembly quality using GiaB HG002

TODO

### Variant detection

We compare our assembled genomes against the published T2T reference genome T2T-CHM13 [Publication] to find genetic variants.  For correct variant calling on sex chromosomes we first determined the sex of the sample using a comparison of UL nanopore read coverages on the T2T-CHM13 X and Y chromosomes. Reads were mapped with _minimap2_[Publication], coverages were calculated using _mosdepth_ [Publication]. A X/Y coverage ratio >0.5 results in male else female assignments. Fully phased small variants and indels were called with _dipcall_ [Publication], larger structural variants with hapdiff[Publication]. We used _whatshap stats_[Publication] to generate variant statistics and estimated phasing accuracy using _whatshap compare_ on a subset of variants appearing present in our samples and the genome in a bottle reference set (CHM13v2.0_HG2-T2TQ100-V1.1) [Publication].

Comment Stephan: Detection of variants compared to hg38, SNVs, indels, SVs, MEIs. Caspar: What are MEIs?

### Haplotype-phased 3D genome using Pore-C

Contact maps were generated using the he _wf-pore-c_ pipeline [Publication]. Individual flowcells were first analysed separately, allowing for granular sequencing quality control. Multiple pairs file from the same sample were merged with _pairtools merge_[Publication]. For haplotype separation of the Pore-C reads we adapted the haplotagging approach used by [Chien et al., 2025] to use the phased SNPs obtained by variant calling from the de novo assembly against the reference genomes. As described in the study we used the Falign alignment method to map the Pore-C reads against the reference genome. We reduced the SNP variant set to high quality regions wnere the mapped Pore-C reads have a minimum depth of 5.  Every chromosome was processed individually, inter-chromosomal contacts were not included in the phasing. We haplotagged the Pore-C reads with _dip3d haplotag_ [Publication] and calculated stats with _dip3d stats_. Finally we merged the individual chromosomes for each sample and extracted separate alignments for every haplotype. For downstream analysis we converted the BAM files to PAF format using _samtools_ and then used a custom script to convert PAF to Pairs format. Now we were able to use the same downstream analysis routines for both the haplotype separated and merged contact pairs files.

poreC, TAD computation, phasing, differential TAD analysis between haplotypes, QTL - variants affecting TADs between haplotypes

### Haplotype-phased deep Methylomes

Combine all libraries (normal and ultra-long) to generate deep, haplotype-phased reference methylomes

### Population structure / ancestry analysis

Compare variants to HapMap / 1000GP, PLINK …
Amounts of variants in various populations

### Pangenome graph

Pan Genome graph.

### Centromere and telomere structures

Analysis with various tools of the repeat composition, length etc.
Repeat composition of each telomere, comparison between the different telomeres, and differences between individuals

### Dark genome - duplicate genes, HLA etc.

We should include one example of a duplicate gene cluster/difficult region that can be used with. In T2T-CHM13 for example they looked a the FRG1 paralogs. 

## Results

### Identification of T2T genomes

T2T assembly quality was determined by integrating telomeric motif detection, chromosome alignment, and gap annotation data for each sample. Contigs were classified as T2T candidates if they contained telomeric repeats at both ends and successfully aligned to reference chromosomes. Assembly quality categories were defined based on remaining gaps: Complete (T2T, zero gaps), Single Small/Large Gap (one gap ≤/>50 kb), Multiple Small/Large Gaps (multiple gaps totaling ≤/>100 kb), or Not Found (no T2T contig identified). Gap statistics were extracted from BED format files and cross-validated against alignment data.



TODO
Structure see above

### Nanopore-only assemblies

- Key QC metrics:
  - Number of perfect chromosomes
  - Number of T2T chromosomes with small issues
  - Length
  - Completeness
- Minimum number of flowcells (Supp Figure)
- Optimal assembly approach (Supp Figure)

## Discussion

Data will be made availabe, can be used as background data or for more detailed analyses.

## References

Tools:
- Mosdepth
- Minimap2
- Dipcall
- Whatshap
- wf-pore-c

GIAB reference set CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz


## Acknowledgement
This work was supported by grants via the BEGIN program of Baden-Württemberg to O.R. and by the European Commission and the BMFTR via the “Genomes of Europe” initiative to O.R.. We thank Ms. Loitz for assistance in blood sampling and …….

## Tables

## Figures

### Figure 1: Assembly quality

a) For every chromosome, Number of completely assembled haplotypes. Barplot.
b) Assembly QC Plot for a single sample showing all gaps and assembly errors
c) Gap Profiling (Location, length), Assembly Error profiling (Classification, Location)
d) Consistency of assembly QC. Correlation Read quality, amount to Assembly  QC, maybe include the subsampling approach here?

### Figure 2: Population genetics

a) SNP PCA (with 1000G background)
b) Worldmap with participants origin
c) Pangenome graph

### Figure 3: 3D Genome

a) TADs heatmap of X chr region, phased, unphased with loops and methylation.
b) Number of contacts for all samples
c) Try to find 3D Structure (phased) with ASHIC

#### Figure 4: Centromere

a) Phased centromer alpha sat plot.
b) Variability of 

## Supplements

### Table 1/Figure 1:

Optimal number of flowcells, subsampling approach. 2 is enough, where does it get better?

### Table 2: Assembly method comparison

QC Results for 1 sample, comparing Verkko (Default), Hifiasm ONT, Hifiasm Herro, Verkko + Hifiasm with AssemblyRepairer

### Table 1: Full assembly QC

Full table of all assembly QC Metrics, all 20 samples.

### Table 2: PoreC quality control

MultiQC Report for individual flowcells.
Merged table with QC per sample.



### Table 3: Correctly assembled centromeres

Results from 
