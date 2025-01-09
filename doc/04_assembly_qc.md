# Assembly QC
Documentation of used QC parameters. Comparison of different assembly and scaffolding methods.

## Collapsed misassemblies

Implemented in this pipeline: paftools asmgene. 
Output described in a [blog post](https://lh3.github.io/2020/12/25/evaluating-assembly-quality-with-asmgene) by Heng Li

MMC = Missing multi copy genes

## Pipeline description Nanopore ONT

![Nanopore T2T Pipeline](img/assm-workflow-1.png)


## Strategies described for assembly polishing:

https://github.com/arangrhie/T2T-Polish/tree/master


## Source 

How to use GFAse with Verkko:

Homopolymer compressed graph needs to be decompressed. 
- `12_uncompress.py`: [GitHub](https://raw.githubusercontent.com/skoren/verkkohic/refs/heads/master/uncompress.sh) S. Koren ([Discussion](https://github.com/marbl/verkko/issues/302#issuecomment-2485843127)) 

We still need to compare the performance of GFAse scaffolding compared to verkko only. 

## Assembly Metrics

### Missassemblies/Collapsed regions
asmgene

### 

### Base quality
Questions do we want to sequence one reference sample (HG002?) to estimate the base error rate?

> Base quality was estimated using yak41, based on the k-mer content of Illumina short reads. Each phased assembly was evaluated separately. K-mer in the short reads were counted using “yak count -b 37”, and quality values (QV) were estimated using “yak qv -K 3.2g -l 100k”. For HG002, we used the 30x Illumina Novaseq PCR-free read set publically available at the Google bucket gs://deepvariant/benchmarking/fastq/wgs_pcr_free/30x/. For the 4 samples from the HPRC (HG01993, HG02132, HG02647, and HG03669), we used 30x Illumina short-reads from the high coverage readset of the 1000 Genomes Project samples (Lorig-Roach et al. 2023)

### Identify T2T contigs

We use the https://github.com/prasad693/Tel_Sequences repo to identify complete T2T chromosomes. 

Steps: 
1) Extract sequence information with seqkit. File: '$output.seqinfo.txt'. 3 Column tsv with: `ID \t seq_length \t 0`
2) Take 1000 bases from start and end with `seqkit subseq` and write into `$output.teloinput.fa`. Append _Start and _End to fasta header.
3) Search telmere motif with tidk. Output file: '$output.haplotype1_telomeric_repeat_windows.csv' Columns: `id,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat`
4) Reformat tidk output to `$output_motif_T2T.txt`. It's filtering for regions with sufficient telomere occurences (≥15 in columns 3 and 4), extracting contig information, and identifying sequences that appear exactly twice in the dataset (found on both ends). Columns: `id \t len`
5) Use minigraph to map ref to assembly and use cov_cal -T to detect T2T contigs. Not additional where this tools comes from. Output: `$output_alignment.txt` 
6) Convert output to new file with the follwing format: 