# 3D Genome

Collection of additional informations related to 3D genome, Dip3D, HiCExplorer.
Need good quality female sample.

## Relation het rate to Pore-C phasing 

Example: TUE_T2T00

      Variants in VCF:    5094598
         Heterozygous:    3562995    (   2998196    SNVs)
               Phased:    3503359    (   2944048    SNVs)
  Heterozygous phased:         98.3% (        98.2% SNVs)
    Heterozygous rate:         70%
  Heterozygosity:    1.1 SNV / 1kb

In Chen et al. 2025:

- Mouse model 7 SNV/1kb
- HG001 <1 SNV/1kb

## Haplotype differences in 3D Genome

Chen et al.:

- Entire X chromosome, X-chromosome inactivation (XCI)-related allele specific expression (ASE)
- X Chr Superdomain at ~115Mb with supterloop formed between tandem repeats DXZ4 and FIRRE, this region also shown in ASHIC. This region also in O-Pore-C. (Deshpande et al., Identifying synergistic high-rder 3D fro nanopore conatemer, Nat Biotechnol. 2022) Froberg et al 2018, Nature comm.
- HIDAD (chr11:1.50 - 2.40Mb): allele-specific subTAD organization at Dlk1-Dio3 and Igf2-H19 (Lleres et al 2019)
- chrX 42Mbp - 56Mbp


## Metrics and analyses currently missing

- SDPR = (Count of distal fragments / Count of proximal fragment). Measure of chromosome compartion. Needs monomer groups with no fewwer then five monomers.
- Monomer count per read

Chen et al used 190x Pore-C reads and were unable to chromsome interaction based imprinting effetcs on autosomes. Check if we can find something for this, else focus on XCI. They only looked at mice and Celllines.

## Brainstorming additional PoreC Analysis

These analyses probablt dont fit in the the first paper.

- Use different phasing approach, replace dip3d. Idea: Map Monomers to the diploid genome. If the reads map uniquely they should be already phased. Use or adapt the haplotype imputation method from Dip3D (3 Steps). 
- Is there any information in the trans contacts? 
- Extended QC, QuasarQC 
