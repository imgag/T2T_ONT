# T2T ONT

Project repository

## Part1: Test pipeline with published datasets

From Koren et al 2024:
- UL:  https://s3-us-west-2.amazonaws.com/humanpangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/nanopore/ultra-long/
- Duplex: https://s3-us-west2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/
- Pore-C: https://www.ncbi.nlm.nih.gov/sra/?term=SRR27664048
- Simplex: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/simplex/

## T2T Reference genome

Used for mapping and comparison of final assemblies: [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)

Annotation for checking genome completeness:
https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/
Ensembl cdna dataset. The header is for GRCh38 alignment, i havent found T2T based annotation yet. cDNA should still be valid, so can be used with T2T after new alignment.
There is also a [comparison](https://ccb.jhu.edu/T2T.shtml) of the RefSeq alignment to both GRCh38 and T2T.

## Assembly polishing 

There is a 3rd method to generate HQ assembly reads:
[Nanopore Polishing Kit](https://labs.epi2me.io/lc2024_t2t/)


## SUP basecalling


### 24070

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
running

#### 04400
```
bin/dorado-0.8.3-linux-x64/bin/dorado basecaller \
    data/ref/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    /mnt/storage3/raw_data/VULCAN/24070/24070LRa010_04400/20241111_1433_2B_PBA36199_3fb5db31 \
    --recursive \
    --trim all \
    --output-dir data/basecalled/24070LRa010_04400 
```

todo