# Methods

Analysis methods. Collection of informations in instructions for tools used in the pipeline. Comparison and evualuation of different analysis approaches. Links to interesting methods that could be used in the project. Small code snippets.

## Use Herro corrected reads for assembly

[Reply](https://github.com/marbl/verkko/issues/322) from S.Koren about best practive for Herro input to T2T

> You can see some info on the README as well as #321. We usually use longest 50x or whatever can be corrected in a reasonable time with Herro. The suggestions pertaining to minimum k-mer coverage in #321 can be used to address the increasing complexity in the graph with high coverage. There are uncorrected regions with systematic errors that start confirming each other and need to be filtered. You can also try the hifiasm correction as an alternative to herro.
> You should not input corrected reads as UL data. The correction can trim and lose parts of the sequence so you want to allow verkko to fill those gaps with the original UL data. Read length is fine, longer reads will build a more resolved initial graph.

## Strategies described for assembly polishing

<https://github.com/arangrhie/T2T-Polish/tree/master>

### Scaffolding

<https://github.com/zengxiaofei/HapHiC>
YaHS

## RepeatMasker

Download DFAM database and put into the environment folder `/mnt/storage2/users/ahgrosc1/environments/envs/snakemake/envs/repeatmasker/share/RepeatMasker/Libraries/famdb$`

Configure RepeatMasker (Use HMMER3.1 & DFAM as search engine)

## RFHAP

Manual Execution on TUE02_03UL dataset, followed by haploid assembly with hifiasm (primary)
Very good separation of haplotypes, also with Herro corrected reads

### 02_3UL: Child reads UL

This is the correct folder. For separation of HERRO corrected reads into the two haplotyped sets.

### 03_3UL_mod: Separate UL bam with methylation into haplotypes

Use assembled haplotypes (.fasta) as input for kmer database.
Seems to have worked

### 04_3UL_trio: Assign correct haplotypes to assembled genomes

Idea:

- Create kmer database from haplotypes
- A) Assign origin for assembled contigs (should be fast)
- B) Generate chunks of assembled contigs and assign haplotype individually

### Hifiasm

Compare following options:

1) *hifiasm_rfhap*: Run hifiasm individually on the separated haplotype reads
2) *hifiasm_ont*: Include poreC scaffolding for the hifiasm
3) *hifiasm_trio*: Use trio binning (kmers) integrated into hifiasm
4) *hifiasm_trio_herro*: Use trio
