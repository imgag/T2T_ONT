# Tasks

## Basecalling

- [ ] Redo 25006LRa108_T2T14_05592 (+Backup)
- [ ] PoreC batch mid August

## Assembly

- [ ] Benchmark Hifiasm with YAHS scaffolding,= use synthetic pairs output from _wf-pore-c_
- [ ] Speed Up Verkko Scaffolding step.
- [ ] Assembly T2T0

## Assembly Quality

- [ ] Redo hifiasm QC table with latest version of scripts
- [ ] Make downsampling plots and stats for all new QC, including Gaps. Create 1 summary figure/table for talks or publication.
- [ ] Add heterozygosity to QC table
- [ ] Add summary_stats QC for SUP basecalled reads / MultiQC. Aggregate for every sample.
- [ ] Compare our samples to HPGRC pipeline
- [ ] Check the accrocentric chromosomes, are there any that were assembled gap free?

## General

- [ ] Add polishing with dorado
- [ ] Check Polishing results
- [ ] Dont calculate mod bases for PoreC!
- [ ] Check where the sequences not in haplotype1 and haplotype2 are coming from.
- [ ] Create liftover ASM -> T2T-CHM13
- [ ] Cleanup analysis_other folder, add required to Git, delete not needed
- [ ] Decide what to do with rfhap, is it required? Maybe for methylation?
- [ ] Create Report for every assembly. Maybe Typst? Markdown? Key metrics and Figures for every sample.

## Base modifications

- [ ] Better implementation of the methylation pipeline.
- [ ] Map UL to T2T-CHM13, pileup
- [ ] Map UL to ASM, pileup
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

## Figures

- [ ]

## Writing

- [ ] Add documentation header to all scripts
- [ ] Remove unused scripts
