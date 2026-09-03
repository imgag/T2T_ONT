# Assembly QC


> "For an assembly to be truly T2T, it must both cover the whole of each chromosome without gaps and be free from large-scale assembly errors. It is critical to rigorously assess the quality of the assembly before concluding it is T2T."
>
> — Li and Durbin, 2024. *Genome Assembly in the Telomere-to-Telomere Era.* Nature Reviews Genetics

Documentation of QC metrics used to evaluate assembly quality.

## QC values overview

### Paftools.js

1. **Length**: Total length of the assembly.
2. **NG50**: Length of the shortest contig at 50% of the total genome length.
3. **NG75**: Length of the shortest contig at 75% of the total genome length.
4. **NGA50**: Length of the shortest aligned contig at 50% of the total genome length.
5. **Rcov**: Percentage of reference genome (T2T-CHM13) covered by the assembly.
6. **Qcov**: Percentage of the query (de novo assembly) covered by the reference.
7. **Rdup**: Percentage of reference genome duplicated in the assembly.
8. **MMC**: Missing Multi-Copy genes fraction, calculated as 1 - |{MCinASM} ∩ {MCinREF}| / |{MCinREF}|.
9. **genome_completeness**: Percentage of single-copy genes present in the assembly compared to the reference.

### Base error

1. **error_rate**: Error rate in the assembly, calculated with k-mers in Merqury.
2. **qv**: Quality value score. -10 × log10(error_rate).

### Whatshap phasing statistics

1. **block_n50**: N50 of phased blocks.
2. **blocks**: Total number of phased blocks.
3. **bp_per_block_avg**: Average number of base pairs per phased block.
4. **bp_per_block_max**: Maximum base pairs in a single phased block.
5. **bp_per_block_median**: Median base pairs per phased block.
6. **bp_per_block_min**: Minimum base pairs in a phased block.
7. **bp_per_block_sum**: Total base pairs across all phased blocks.
8. **heterozygous_snvs**: Number of heterozygous SNVs.
9. **heterozygous_variants**: Total number of heterozygous variants.
10. **phased**: Number of phased variants.
11. **phased_fraction**: Fraction of heterozygous variants that have been phased.
12. **phased_snvs**: Number of heterozygous SNVs that have been phased.
13. **phased_snvs_fraction**: Fraction of heterozygous SNVs that have been phased.
14. **singletons**: Number of singleton variants.
15. **unphased**: Number of heterozygous variants not phased.
16. **variant_per_block_avg**: Average number of variants per phased block.
17. **variant_per_block_max**: Maximum number of variants in a single phased block.
18. **variant_per_block_median**: Median number of variants per phased block.
19. **variant_per_block_min**: Minimum number of variants in a phased block.
20. **variant_per_block_sum**: Total number of variants across all phased blocks.
21. **variants**: Number of biallelic variants in the input VCF excluding duplicate positions.

## Detailed metric descriptions

### Collapsed misassemblies (MMC/MSC)

Implemented via `paftools asmgene`. Described in a [blog post](https://lh3.github.io/2020/12/25/evaluating-assembly-quality-with-asmgene) by Heng Li.

- **MMC** = Missing Multi-Copy genes
- **MSC** = Missing Single-Copy genes

> "paftools asmgene" detects missing genes by aligning transcripts to both an assembly haplotype and a haploid reference and counting discrepancies in gene copy number. Subsequently, the percentage of genes that are multi-copy in the haploid reference but not in the assembly haplotype (%MMC) and the percentage of genes that are single-copy in the haploid reference but not in the assembly haplotype (%MSC) were computed.

### T2T contig identification

Complete telomere-to-telomere contigs are identified using the [Tel_Sequences](https://github.com/prasad693/Tel_Sequences) approach:

1. Extract sequence information with seqkit: `seqkit fx2tab -j $threads -C N -l -n -o $output.seqinfo.txt $assembly`
2. Extract 1000 bases from start and end with `seqkit subseq`.
3. Search for telomere motifs with `tidk`. Sequences with ≥15 telomeric repeats on both ends are classified as T2T candidates.
4. Map reference to assembly with `minigraph` and use `cov_cal -T` to confirm ≥95% coverage on both query and reference in a single contiguous alignment.
