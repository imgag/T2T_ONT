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

### problem with hic for diploid
- why was it never a problem with haploid? 
1.  no trans contacts
2. all contacts are sorted by chr already bc the bam has been merged from phased single chr bams = pairs already in chromosomal order (1-1,2-2 etc)

3. maybe this can be fixed easier by doing pairtools parse2 --flip and paitools sort

- 25006LRa064_T2T06_PoreC_05350.ns.bam name sorted bam for each fc
- format: 
893dd98b-36a7-41b5-8b91-0e5aa27e5578:0000:0329
893dd98b-36a7-41b5-8b91-0e5aa27e5578:0329:0336
893dd98b-36a7-41b5-8b91-0e5aa27e5578:0336:0864
.....

- todo: count frags per read, count frag size, mapped frags per read

-  this might NOT be correct!!! the ns bam is a paired-end. need to check 1 read. - but this checks out??
-  why new nonadj bam so much larger than old (wf-pore-c also use expand? is it bc of restriction filters?)

rule flip_pairs_wfporec:
    input:
        pairs = "../analysis_other/wf-pore-c/25006LRa064_T2T06_PoreC_05350/pairs/25006LRa064_T2T06_PoreC_05350.nonadj.pairs.gz",
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        flipped = "outputs/wf-pore-c/25006LRa064_T2T06_PoreC_05350/pairs/25006LRa064_T2T06_PoreC_05350.nonadj_flipped.pairs.gz"
    conda:
        "../env/pairtools.yml"
    log:
        "logs/flip_pairs_25006LRa064_T2T06_PoreC_05350.log"
    shell:
        """
        pairtools flip \
        --chroms-path {input.chromsize} \
        --output {output.flipped} {input.pairs}
        """

    
rule clean_pairs2:
    input:
        "outputs/wf-pore-c/25006LRa064_T2T06_PoreC_05350/pairs/25006LRa064_T2T06_PoreC_05350.nonadj_flipped.pairs.gz"
    output:
        temp("outputs/hic_files_T2T/25006LRa064_T2T06_PoreC_05350/pairs_juice/25006LRa064_T2T06_PoreC_05350_flipped.pairs.for_juice")
    shell:
        """
        zcat {input} | \
        grep -v '^#' | \
        awk '{{OFS="\\t"; print $1,$6,$2,$3,0,$7,$4,$5,1,$11,$12}}' | \
        awk '{{OFS="\\t"; $2=($2=="+")?0:1; $6=($6=="+")?0:1; print}}' | \
        sort -k3,3V -k7,7V > {output}
        """
- for parse2:
    - Until recently, I'd been using pairtools parse. Recently I switched to parse2. In testing my parse2-based pipeline, I started running into issues with sorting where, for example, I'd have a (chr1, chr16) pair in one section of the .pairs file and (chr16, chr1) in another section. The problem appears to be resolved now that I started using the --flip parameter.
    - I don't think this was ever an issue when I used parse, and I don't think I used the --flip parameter with parse. So I'm wondering if possibly pairtools has a different default value for the --flip parameter for parse vs. parse2? The documentation isn't clear on this (and could stand to be a bit more clearly worded, as this is a bit ambiguous:

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

