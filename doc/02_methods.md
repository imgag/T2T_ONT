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
Seems to have worked

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


```
Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>
Options:
  Input/Output:
    -o STR       prefix of output files [hifiasm.asm]
    -t INT       number of threads [1]
    -h           show help information
    --version    show version number
  Preset options:
    --ont        assemble Oxford Nanopore reads
  Overlap/Error correction:
    -k INT       k-mer length (must be <64) [51]
    -w INT       minimizer window size [51]
    -f INT       number of bits for bloom filter; 0 to disable [37]
    -D FLOAT     drop k-mers occurring >FLOAT*coverage times [5.0]
    -N INT       consider up to max(-D*coverage,-N) overlaps for each oriented read [100]
    -r INT       round of correction [3]
    -z INT       length of adapters that should be removed [0]
    --max-kocc   INT
                 employ k-mers occurring <INT times to rescue repetitive overlaps [2000]
    --hg-size    INT(k, m or g)
                 estimated haploid genome size used for inferring read coverage [auto]
  Assembly:
    -a INT       round of assembly cleaning [4]
    -m INT       pop bubbles of <INT in size in contig graphs [10000000]
    -p INT       pop bubbles of <INT in size in unitig graphs [0]
    -n INT       remove tip unitigs composed of <=INT reads [3]
    -x FLOAT     max overlap drop ratio [0.8]
    -y FLOAT     min overlap drop ratio [0.2]
    -i           ignore saved read correction and overlaps
    -u           post-join step for contigs which may improve N50; 0 to disable; 1 to enable
                 [1] and [1] in default for the UL+HiFi assembly and the HiFi assembly, respectively
    --hom-cov    INT
                 homozygous read coverage [auto]
    --lowQ       INT
                 output contig regions with >=INT% inconsistency in BED format; 0 to disable [70]
    --b-cov      INT
                 break contigs at positions with <INT-fold coverage; work with '--m-rate'; 0 to disable [0]
    --h-cov      INT
                 break contigs at positions with >INT-fold coverage; work with '--m-rate'; -1 to disable [-1]
    --m-rate     FLOAT
                 break contigs at positions with <=FLOAT*coverage exact overlaps;
                 only work with '--b-cov' or '--h-cov'[0.75]
    --primary    output a primary assembly and an alternate assembly
    --ctg-n      INT
                 remove tip contigs composed of <=INT reads [3]
  Trio-partition:
    -1 FILE      hap1/paternal k-mer dump generated by "yak count" []
    -2 FILE      hap2/maternal k-mer dump generated by "yak count" []
    -3 FILE      list of hap1/paternal read names []
    -4 FILE      list of hap2/maternal read names []
    -c INT       lower bound of the binned k-mer's frequency [2]
    -d INT       upper bound of the binned k-mer's frequency [5]
    --t-occ      INT
                 forcedly remove unitigs with >INT unexpected haplotype-specific reads;
                 ignore graph topology; [60]
    --trio-dual  utilize homology information to correct trio phasing errors
  Purge-dups:
    -l INT       purge level. 0: no purging; 1: light; 2/3: aggressive [0 for trio; 3 for unzip]
    -s FLOAT     similarity threshold for duplicate haplotigs in read-level [0.75 for -l1/-l2, 0.55 for -l3]
    -O INT       min number of overlapped reads for duplicate haplotigs [1]
    --purge-max  INT
                 coverage upper bound of Purge-dups [auto]
    --n-hap      INT
                 number of haplotypes [2]
  Hi-C-partition:
    --h1 FILEs   file names of Hi-C R1  [r1_1.fq,r1_2.fq,...]
    --h2 FILEs   file names of Hi-C R2  [r2_1.fq,r2_2.fq,...]
    --seed INT   RNG seed [11]
    --s-base     FLOAT
                 similarity threshold for homology detection in base-level;
                 -1 to disable [0.5]; -s for read-level (see <Purge-dups>)
    --n-weight   INT
                 rounds of reweighting Hi-C links [3]
    --n-perturb  INT
                 rounds of perturbation [10000]
    --f-perturb  FLOAT
                 fraction to flip for perturbation [0.1]
    --l-msjoin   INT
                 detect misjoined unitigs of >=INT in size; 0 to disable [500000]
  Ultra-Long-integration:
    --ul FILEs   file names of Ultra-Long reads [r1.fq,r2.fq,...]
    --ul-rate    FLOAT
                 error rate of Ultra-Long reads [0.2]
    --ul-tip     INT
                 remove tip unitigs composed of <=INT reads for the UL assembly [6]
    --path-max   FLOAT
                 max path drop ratio [0.6]; higher number may make the assembly cleaner
                 but may lead to more misassemblies
    --path-min   FLOAT
                 min path drop ratio [0.2]; higher number may make the assembly cleaner
                 but may lead to more misassemblies
    --ul-cut     INT
                 filter out <INT UL reads during the UL assembly [0]
  Dual-Scaffolding:
    --dual-scaf  output scaffolding
    --scaf-gap   INT
                 max gap size for scaffolding [3000000]
  Telomere-identification:
    --telo-m     STR
                 telomere motif at 5'-end; CCCTAA for human [NULL]
    --telo-p     INT
                 non-telomeric penalty [1]
    --telo-d     INT
                 max drop [2000]
    --telo-s     INT
                 min score for telomere reads [500]
  ONT simplex assembly (beta):
    --ont        assemble ONT simplex reads in fastq format
    --chem-c     INT
                 detect chimeric reads with <=INT other reads support [1]
    --chem-f     INT
                 length of flanking regions for chimeric read detection [256]
    --rl-cut     INT
                 filter out ONT simplex reads shorter than <INT> for assembly [1000]
    --sc-cut     INT
                 filter out ONT simplex reads with a mean base quality score below <INT> [10]
Example: ./hifiasm -o NA12878.asm -t 32 NA12878.fq.gz