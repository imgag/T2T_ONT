# Methods

Analysis methods. Collection of informations in instructions for tools used in the pipeline. Comparison and evualuation of different analysis approaches. Small code snippets.

## Use Herro corrected reads for assembly

[Reply](https://github.com/marbl/verkko/issues/322) from S.Koren about best practive for Herro input to T2T 

> You can see some info on the README as well as #321. We usually use longest 50x or whatever can be corrected in a reasonable time with Herro. The suggestions pertaining to minimum k-mer coverage in #321 can be used to address the increasing complexity in the graph with high coverage. There are uncorrected regions with systematic errors that start confirming each other and need to be filtered. You can also try the hifiasm correction as an alternative to herro.
> You should not input corrected reads as UL data. The correction can trim and lose parts of the sequence so you want to allow verkko to fill those gaps with the original UL data. Read length is fine, longer reads will build a more resolved initial graph.

## Strategies described for assembly polishing:

https://github.com/arangrhie/T2T-Polish/tree/master


### Basecalling

Basecalling is not yet integrated into the pipeline. Should be done when all samples are assembled in routine.
*Include MOD basecalling for UL reads in next basecalling*

#### #04339
```
bin/dorado-0.8.3-linux-x64/bin/dorado basecaller \
    data/ref/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    data/raw/set1/24070/24070LRa003_04399/20241111_1425_2A_PBA32113_7de81dc2 \
    --recursive \
    --trim all \
    --output-dir data/basecalled/24070LRa003_04399

bin/dorado-0.8.3-linux-x64/bin/dorado summary \
    data/basecalled/24070LRa003_04399/*.bam \
    > data/basecalled/24070LRa003_04399/sequencing_summary.txt
```

#### 04400
```
bin/dorado-0.8.3-linux-x64/bin/dorado basecaller \
    data/ref/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    /mnt/storage3/raw_data/VULCAN/24070/24070LRa010_04400/20241111_1433_2B_PBA36199_3fb5db31 \
    --recursive \
    --trim all \
    --output-dir data/basecalled/24070LRa010_04400 
```
#### Scaffolding

https://github.com/zengxiaofei/HapHiC


#### RepeatMasker

Download DFAM database and put into the environment folder `/mnt/storage2/users/ahgrosc1/environments/envs/snakemake/envs/repeatmasker/share/RepeatMasker/Libraries/famdb$`

Configure RepeatMasker (Use HMMER3.1 & DFAM as search engine):

#### RFHAP

Manual Execution on TUE02_03UL dataset

##### 01_3UL_childherro: Child reads HQ herro corrected

```
(rfhap) ahgrosc1@SRV026:/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT$ bash analysis_other/rfhap/rfhap.sh 

 N E X T F L O W   ~  version 24.10.5

Launching `bin/rfhap-1.0/rfhap.nf` [loving_brahmagupta] DSL2 - revision: 5445b850d9

executor >  local (69)
[aa/b4c9c0] create_kmers_database (4) [100%] 4 of 4 ✔
[2f/7abcc0] print_paths (4)           [100%] 4 of 4 ✔
[a9/9fe7e6] sort_kmer_db (1)          [100%] 1 of 1 ✔
[1e/f085f3] FastKM (1)                [100%] 1 of 1 ✔
[3a/de5571] trainRF (1)               [100%] 1 of 1 ✔
[b9/d7715d] predictRF (28)            [100%] 54 of 54 ✔
[b5/d8a2a0] setHaplotypes (1)         [100%] 1 of 1 ✔
[14/ecea86] seqtk (1)                 [100%] 3 of 3 ✔
Completed at: 24-Apr-2025 09:59:59
Duration    : 12h 25m 44s
CPU hours   : 543.5
Succeeded   : 69
```
Output file sizes: 

```  
 39G Apr 24 09:59 TUE_02_03UL.HQ_herro.fastq_fastkm_matrix.hapA.fq.gz
 38G Apr 24 09:59 TUE_02_03UL.HQ_herro.fastq_fastkm_matrix.hapB.fq.gz
2.9G Apr 24 09:45 TUE_02_03UL.HQ_herro.fastq_fastkm_matrix.hapU.fq.gz
```

Very good separation of haplotypes, also with Herro corrected reads

##### 02_3UL: Child reads UL

As they should be

##### 03_3UL_mod: Separate UL bam with methylation into haplotypes

Use assembled haplotypes (.fasta) as input for kmer database. 
.. running

##### 04_3UL_trio: Assign correct haplotypes to assembled genomes

Idea:
- Create kmer database from haplotypes
- A) Assign origin for assembled contigs (should be fast)
- B) Generate chunks of assembled contigs and assign haplotype individually


#### Hifiasm

Compare following options:

1) *hifiasm_rfhap*: Run hifiasm individually on the separated haplotype reads
2) *hifiasm_porec*: Include poreC scaffolding for the hifiasm
3) *hifiasm_trio*: Use trio binning (kmers) integrated into hifiasm

