# Tasks

## Basecalling

- [X] Redo 25006LRa108_T2T14_05592 (+Backup)
- [X] PoreC batch mid August
- [ ] P02 in diagnostics

## Assembly

- [ ] Benchmark Hifiasm with YAHS scaffolding, use synthetic pairs output from _wf-pore-c_
- [X] Speed Up Verkko Scaffolding step.

## Assembly Quality

- [ ] Redo hifiasm QC table with latest version of scripts
- [ ] Check out [BISER](https://github.com/0xTCG/biser/) tool for finding segmental duplication
- [X] Make downsampling plots and stats for all new QC, including Gaps. Create 1 summary figure/table for talks or publication.
- [ ] Add heterozygosity to QC table
- [ ] Add summary_stats QC for SUP basecalled reads / MultiQC. Aggregate for every sample.
- [X] Compare our samples to HPGRC pipeline
- [X] Check the accrocentric chromosomes, are there any that were assembled gap free?
- [X] Add hmm-flagger QC
- [X] Rename repeatmasker and nucflag outputs
- [X] Add nucflag and repeatmaskter outputs to Git

## General

- [X] Add polishing with dorado
- [X] Select polishing files with dorado
- [X] Check Polishing results
- [X] Dont calculate mod bases for PoreC!
- [X] Check where the sequences not in haplotype1 and haplotype2 are coming from.
- [X] Create liftover ASM -> T2T-CHM13
- [ ] Cleanup analysis_other folder, add required to Git, delete not needed
- [ ] Decide what to do with rfhap, is it required? Maybe for methylation?
- [ ] Create Report for every assembly. Maybe Typst? Markdown? Key metrics and Figures for every sample.
- [ ] Upload folder: Generate README describing files, update file names

## Base modifications

- [X] Better implementation of the methylation pipeline.
- [X] Map UL to T2T-CHM13, pileup
- [X] Map UL to ASM, pileup
- [ ] Check if we want to use rfhap

## 3D Genome

- [ ] Create figure with multiple haplotypes
- [ ] Calculate heterozygosities for every sample (reference to ancestry) and compare to performance of PoreC phasing. Assumption: Higher heterozygosity allows better haplotype separation. Baseline HG001 (low heterozygostiy) <1 het-SNV per 1kb
- [ ] Additional QC (QuasarQC)
- [ ] Check on imprinted alleles
- [ ] Phase reads mapped on assembly instead of reference genome
- [ ] Idea: Adapt haplotype imputation process for PoreC mapping to diploid reference genome. Compare against Dip3D paper.
- [ ] Use ASHIC to compare 3d structure, Chen et al: X-specific bipartite structure at DXZ4

## Centromeres

- [ ] GCP: Read paper in detail, apply to our results
- [ ] Create GQC plots

## Population genetics

- [ ] Ancestry estimation STRUCTURE/ADMIXTURE

## Figures

- [ ]

## Writing

- [ ] Add documentation header to all scripts
- [ ] Describe all final QC Parameters
- [ ] Remove unused scripts
