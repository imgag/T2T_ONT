# Tasks

## Basecalling

- [X] Redo 25006LRa108_T2T14_05592 (+Backup)
- [X] PoreC batch mid August

## Assembly

- [ ] Create PoreC scaffolded Hifiasm assembly, use script from ONT or create new
- [ ] Trio Hifiasm
- [ ] Trio Verkko

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
- [ ] Create Upload folder: Generate README describing files, update file names. Options: IMGAG Cloud or separate S3:// instance. Either on imgag.de or external storage.
- [ ] Create Elrin GitHub Account with data overview and download links


## Base modifications

- [X] Better implementation of the methylation pipeline.
- [X] Map UL to T2T-CHM13, pileup
- [X] Map UL to ASM, pileup
- [ ] Check if we want to use rfhap

## 3D Genome

- [ ] Create QC Figure.
    - [ ] Percentage of phased reads
    - [ ] Monomers per fragment?
    - [ ] IGF19 Region for one sample
    - [ ] XCI for one sample

## Centromeres 

- [ ] Create Figure for Cenmap results

## Population genetics

- [ ] Use lifted dataset from GitHub
- [ ] Fix sample bias in PCA
- [ ] Research how many SNVs should be used

## Figures

- [ ] Create figures
- [ ] Create figure annotations

## Writing

- [ ] Add documentation header to all scripts
- [ ] Describe all final QC Parameters
- [ ] Remove unused scripts


## Ideas for later

- [ ] Check out [BISER](https://github.com/0xTCG/biser/) tool for finding segmental duplication
- [ ] Research how many SNVs should be usedp
- [ ] Calculate heterozygosities for every sample (reference to ancestry) and compare to performance of PoreC phasing. Assumption: Higher heterozygosity allows better haplotype separation. Baseline HG001 (low heterozygostiy) <1 het-SNV per 1kb
- [ ] Add polishing with dorado (HQ reads). Does it improve?
- [ ] Investigate APK Polishing. Does it really get worse? Where are the changes located?
