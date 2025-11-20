# Tasks

## Basecalling

- [X] Redo 25006LRa108_T2T14_05592 (+Backup)
- [X] PoreC batch mid August

## Assembly

- [ ] Create PoreC scaffolded Hifiasm assembly, use script from ONT or create new
- [ ] Trio Hifiasm
- [ ] Trio Verkko

## Assembly Quality

- [ ] Decide whether we want to have a karyotype plot with difficult regions.
- [ ] Create some figure with potential assembly issues
- [ ] Compare Nucflag/Flagger results with HG002 samples

## General

- [ ] Check where the sequences not in haplotype1 and haplotype2 are coming from.
- [ ] Create liftover ASM -> T2T-CHM13
- [ ] Cleanup analysis_other folder, add required to Git, delete not needed
- [ ] Decide what to do with rfhap, is it required? Maybe for methylation?
- [ ] Create Report for every assembly. Maybe Typst? Markdown? Key metrics and Figures for every sample.
- [ ] Create Upload folder: Generate README describing files, update file names. Options: IMGAG Cloud or separate S3:// instance. Either on imgag.de or external storage.
- [ ] Create Elrin GitHub Account with data overview and download links


## Base modifications
 
- [ ] Check that PoreC pipeline has finished for all samples, together with Sarah
- [ ] Check if we need modifications against T2T-CHM13
- [ ] Check of modifications are phased correctly, else do that
- [ ] Maybe we could use rfhap for read phasing?

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
