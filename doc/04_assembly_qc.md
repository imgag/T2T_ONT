# Assembly QC
Documentation of used QC parameters. Comparison of different assembly and scaffolding methods.

## QC Values Overview

### When the sample is not a reference genome:

3. **#breaks**: Number of breaks in the assembly. (_asmstat_)
5. **Length**: Total length of the assembly.
6. **NG50**: Length of the shortest contig at 50% of the total genome length.
7. **NG75**: Length of the shortest contig at 75% of the total genome length.
8. **NGA50**: Length of the shortest aligned contig at 50% of the total genome length.
10. **Rcov**: Percentage of Reference genome (T2T) covered by assembly.
9. **Qcov**: Percentage of Query (De novo Assembly) covered by reference .
10. **l_cov**: ??
11. **Rdup**: Percentage of reference genome duplicates in assembly.
19. **n_y_chrom**: Number of Y chromosome sequences.
1. **MMC**: Missing Multi-Copy genes fraction, calculated as 1 - |{MCinASM} ∩ {MCinREF}| / |{MCinREF}|.
2. **genome_completeness**: Percentage of single-copy genes present in the assembly compared to the reference.

### When sample is HG002:

2. **error_rate**: Error rate in the assembly. Calculated with kmers in _Merqury_
3. **qv**: Quality value score. 10^(-error_rate)

### whatshap
1. **block_n50**: The N50 of phased blocks. This is the length of the shortest phased block at 50% of the total phased block length. It provides an indication of the contiguity of the phased blocks.

2. **blocks**: The total number of phased blocks. A phased block is a contiguous segment of the genome where the phase (the arrangement of alleles) is known.

3. **bp_per_block_avg**: The average number of base pairs per phased block. This metric gives an idea of the typical size of phased blocks.

4. **bp_per_block_max**: The maximum number of base pairs in a single phased block. This indicates the size of the largest phased block.

5. **bp_per_block_median**: The median number of base pairs per phased block. This is the middle value of the block sizes when they are sorted in ascending order.

6. **bp_per_block_min**: The minimum number of base pairs in a single phased block. This indicates the size of the smallest phased block.

7. **bp_per_block_sum**: The total number of base pairs across all phased blocks. This gives the overall length of the genome that has been phased.

8. **heterozygous_snvs**: The number of heterozygous single nucleotide variants (SNVs). These are positions in the genome where there are two different alleles present in the individual.

9. **heterozygous_variants**: The total number of heterozygous variants, including SNVs and other types of variants such as insertions and deletions.

10. **phased**: The number of variants that have been phased. Phasing determines which variants are on the same chromosome and which are on the homologous chromosome.

11. **phased_fraction**: The fraction of heterozygous variants that have been phased. This is calculated as the number of phased variants divided by the total number of heterozygous variants.

12. **phased_snvs**: The number of heterozygous SNVs that have been phased.

13. **phased_snvs_fraction**: The fraction of heterozygous SNVs that have been phased. This is calculated as the number of phased SNVs divided by the total number of heterozygous SNVs.

14. **singletons**: The number of singleton variants. These are variants that appear only once in the dataset.

15. **unphased**: The number of heterozygous variants that have not been phased.

16. **variant_per_block_avg**: The average number of variants per phased block. This gives an idea of the density of variants within phased blocks.

17. **variant_per_block_max**: The maximum number of variants in a single phased block.

18. **variant_per_block_median**: The median number of variants per phased block.

19. **variant_per_block_min**: The minimum number of variants in a single phased block.

20. **variant_per_block_sum**: The total number of variants across all phased blocks.

21. **variants**: Number of biallelic variants in the input VCF excluding duplicate positions 


## QC Value detailed description

Explanation of some key metrics, as used in similar publications. 

### Collapsed misassemblies (MMC)

Implemented in this pipeline: paftools asmgene. 
Output described in a [blog post](https://lh3.github.io/2020/12/25/evaluating-assembly-quality-with-asmgene) by Heng Li

MMC = Missing multi copy genes

## Source 

How to use GFAse with Verkko:

Homopolymer compressed graph needs to be decompressed. 
- `12_uncompress.py`: [GitHub](https://raw.githubusercontent.com/skoren/verkkohic/refs/heads/master/uncompress.sh) S. Koren ([Discussion](https://github.com/marbl/verkko/issues/302#issuecomment-2485843127)) 

We still need to compare the performance of GFAse scaffolding compared to verkko only. 

### Base quality
Questions do we want to sequence one reference sample (HG002?) to estimate the base error rate?

> Base quality was estimated using yak41, based on the k-mer content of Illumina short reads. Each phased assembly was evaluated separately. K-mer in the short reads were counted using “yak count -b 37”, and quality values (QV) were estimated using “yak qv -K 3.2g -l 100k”. For HG002, we used the 30x Illumina Novaseq PCR-free read set publically available at the Google bucket gs://deepvariant/benchmarking/fastq/wgs_pcr_free/30x/. For the 4 samples from the HPRC (HG01993, HG02132, HG02647, and HG03669), we used 30x Illumina short-reads from the high coverage readset of the 1000 Genomes Project samples (Lorig-Roach et al. 2023)

### Merqury

Column headers:



### Identification of complete T2T contigs

We use the https://github.com/prasad693/Tel_Sequences repo to identify complete T2T chromosomes. 

Steps: 
1) Extract sequence information with seqkit. File: '$output.seqinfo.txt'. 3 Column tsv with: `ID \t seq_length \t 0`
2) Take 1000 bases from start and end with `seqkit subseq` and write into `$output.teloinput.fa`. Append _Start and _End to fasta header.
3) Search telmere motif with tidk. Output file: '$output.haplotype1_telomeric_repeat_windows.csv' Columns: `id,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat`
4) Reformat tidk output to `$output_motif_T2T.txt`. It's filtering for regions with sufficient telomere occurences (≥15 in columns 3 and 4), extracting contig information, and identifying sequences that appear exactly twice in the dataset (found on both ends). Columns: `id \t len`
5) Use minigraph to map ref to assembly and use cov_cal -T to detect T2T contigs. Not additional where this tools comes from. Output: `$output_alignment.txt` 
6) Convert output to new file with the follwing format: 


### Comparison of flowcell impact

