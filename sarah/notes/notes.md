## ok plan:

### 1. corrected cool of adj + nonadj at suitable resolutions 
--> done

#### 1.1 check .hic, compare to other human hic
--> done
- non adj look better at 25 kb, 10 kb eh, all higher res --> not feasible
- T2T04: no chrx no chry, ? chromsize file? sample specific? old pairs of T2T03 had both.
- roughly 100 million mapped / valid junction reads, is sufficient to support a 40kb data resolution (Lajoie 2016)
---
- nonadj: 337,790,263 (max expansion depth=20)
- adj: 81,349,419 

e.g. d79b2e5f-ac5e-4b64-b6a6-4a47d3f32e5f_0001121076
- 22 non adj pairs
- 7 adj pairs (8 frags?)
---

#### 1.2 add + check plot, decide z-score thresh
--> done:
pairs_to_cooler, merge_mcools, hic_diagnostic_plot, hic_correct_matrix

- #### 30.10.25: changed to cooler balance to stay within cooler, better for downstream cooltools analysis; new: cooler_balance

- how to apply internal threshold calc of tool hicCorrectMatrix?
- way to extract threshold out of diagnostic plot and directly insert into correction function:

hicCorrectMatrix diagnostic_plot -m matrix.h5 -o QC/matrix.pdf &> QC/matrix_mad_threshold.out
madscore=$(grep "mad threshold" QC/matrix_mad_threshold.out | sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g');
upper=$(echo -3*$madscore | bc);
echo $madscore " " $upper > QC/matrix_mad_threshold.values
thresholds=$(cat QC/matrix_mad_threshold.values);
hicCorrectMatrix correct --filterThreshold $thresholds -m matrix.h5 -o matrix_cor.h5

mine:
madscore=$(grep "mad threshold" {input.diagnostic} | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g')
        upper=$(echo -3*$madscore | bc)
        thresholds=$(echo $madscore " " $upper)

find:
    - common number of its for human
    - is ice applied per chrom per default? yes
    - should i normalise the matrices first?

### 2. matrix analysis

A linear correlation between the Pore-C and Hi-C contact 
matrices measured three different similarity metrics:
(1) raw contact matrices, 
(2) compartmental eigenvector scores identified using 
cooltools call-compartments and 
(3) TAD insulation scores calculated using the cooltools 
diamond-insulation tools.

#### 2.1 hic_plot_dist_vs_counts
- compare adj + nonadj, see reduction in sr similar to papers?
- done but idk if it plots relative counts  
--> do own script

#### 2.2 hic_find_tads and hic_detect_loops
- resolution humans? size humans? number?
- does any of them output insulation score? no, but cooltools insulation does
- cooltools diamond-insulation (Desphande 2022)

#### 2.3 NEW:   
- hicCorrelate is a dedicated Quality Control tool that allows the correlation of multiple Hi-C matrices at once with either a heatmap or scatter plots output.
    - --> better used on un-corrected matrices
    - --> --outFileNameScatter

- i will need to select a resolution at which the adj matrix is not shit (results: 40 kb, pcc=0.92, 25 kb pcc=0.88)
- HiCrep v1.2.0: stratum-adjusted correlation coefficient of the pairwise contact matrix between samples (https://github.com/TaoYang-dev/hicrep)

2.4 NEW:
- resolution? --> hg38 at 100kb (cooltools doc)
- In humans and mice, GC content is useful for phasing because it typically has a strong correlation at the 100kb-1Mb bin level with the eigenvector. (this looks very similar to lieberman)
- The cooltools genome command group can generate GC or gene coverage tracks.

- cooltools eigs-cis: very high correlation between adj and non adj at 100 kb
- to do: plot along heatmap. --> done
- done: look at .hic eig --> similar to lieberman 2009 eigs of humans!

 - dCHiC https://github.com/ay-lab/dcHiC: tool for some differential hi-c analysis

##### several forms of balance:
- hicCorrectMatrix, for all downstream hicExplorer tools
- cooler balance, for all downstream cooler and cooltools
- juicer tools pre: hic files with KR can be used for several kinds of tools, e.g. fanc

### 3. general 
- wf-pore-c
    - no expand AND also pairtools restrict, meaning this is much more stringent than the single haps
- maybe this will not be liked by reviewer, aka the "more difficult" datasets are generated with less stringent settings
- best practice may be to treat both same starting from aligned bam outputs


### caspars snakemake
- finished_samples = list of all sample names without hap, can be used in expand

# XCI
- strong insulation changes across the superdomain boundary
- 


# errors#
# juicer tools pre
Hi Marine,

I downloaded your file. Please see below:

awk '{if ($3 == "chr1" && $7 == "chr10") {print NR,$0; exit}}' dijhic006_merge_hg38.input.sorted.txt
4691370 A00417:65:HMMNNDMXX:2:1101:10104:6715    0    chr1    24928772    0    0    chr10    75633923    1    60    47
awk '{if ($3 == "chr10" && $7 == "chr1") {print NR,$0; exit}}' dijhic006_merge_hg38.input.sorted.txt
12783104 A00417:65:HMMNNDMXX:2:1101:10004:2284    16    chr10    70907616    0    16    chr1    12404324    1    60    60

So you have records where chr1 comes before, say, chr10 and those where it is the other way around. This never happens if you use juicer pipeline but can happen with other pipeline. So you first need to make sure that the order of chromosomes for each pair is the same. I would do the following:

awk '{if ($3 <= $7) print $0; else print $1,$6,$7,$8,$9,$2,$3,$4,$5,$11,$10}' dijhic006_merge_hg38.input.sorted.txt > correct_dijhic006_merge_hg38.input.txt
sort -k3,3d -k7,7d correct_dijhic006_merge_hg38.input.txt > correct_dijhic006_merge_hg38.input.sorted.txt
and then run pre on correct_dijhic006_merge_hg38.input.sorted.txt

Hope this helps,
Moshe.

