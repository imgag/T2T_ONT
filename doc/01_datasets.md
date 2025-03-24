# Datasets used in this study

## Trio datasets

Father: GE-MED-T2T07
Mother: GE-MED-T2T08
Child: 2410869

## Published sequencing data

In the first part of the analysis we setup the pipeline and did some initial testing using the published sequencing data. Although we want to strictly generate de novo assemblies without help of a reference genome we are using reference genomes to assess the quality of the assembles genomes. 

From Koren et al 2024:
- UL:  https://s3-us-west-2.amazonaws.com/humanpangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/nanopore/ultra-long/
- Duplex: https://s3-us-west2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/
- Pore-C: https://www.ncbi.nlm.nih.gov/sra/?term=SRR27664048
- Simplex: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/0CB931D5-AE0C-4187-8BD8-B3A9C9BFDADE--UCSC_HG002_R1041_Duplex_Dorado/Dorado_v0.1.1/simplex/

## Reference genomes

Multiple groups are working on a T2T reference genome:
- T2T consortium
- Human Pangenome Reference Consortium (HPRC)
- Genome In A Bottle (GIAB), NIST

Aim: Q100 (1 error in 10 billion bases)

### Haploid reference genome (CHM13-T2T)

The latest release (T2T-CHMv2.0) includes a Y chromosome. Used for mapping and comparison of final assemblies: [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)

Haploid cellline from a hydatidiform mole. This occurs when a sperm fertilizes an egg that lost its DNA. After mitosis the genotype is 46,XX but the genome consists only of a duplicated set of haploid chromosomes.

### Diploid reference genome (HG002)

Published reference genome on [Github](https://github.com/marbl/HG002?tab=readme-ov-file)

Current version is hg002v1.1 (July, 2024)
[Download](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz)

Methods used to polish and correct errors are described in this paper:
Mc Cartney et al. Nature Methods (June 2022): 687–95. https://doi.org/10.1038/s41592-022-01440-3.

There is a draft version of a GIAB benchmarkset using T2TQ100 (HG002) called against CHM13v2.0 here:
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/CHM13v2.0_HG2-T2TQ100-V1.1.vcf.gz

We use this file to analyse our phasing performance.

### Genome annotation

Annotation for checking genome completeness:
https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/
Ensembl cdna dataset. The header is for GRCh38 alignment, i havent found T2T based annotation yet. cDNA should still be valid, so can be used with T2T after new alignment.
There is also a [comparison](https://ccb.jhu.edu/T2T.shtml) of the RefSeq alignment to both GRCh38 and T2T.

## Additional datasets:

ONT Official Blogpost: 
 45X UL -> Herro
 35X PoreC 
 35X 6B4 Polishing



 