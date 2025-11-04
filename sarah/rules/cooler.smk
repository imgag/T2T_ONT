datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

#  hic map analysis with cooler and cooltools 
#  uses cooler balance (ice) for matrix correction as cooltools is used for downstream analyses

rule collect_cooler:
    input:
        expand("outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}.mcool",
        dataset=datasets,
        exp=config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/qc/correlation_{plot}_{resolution}.png",
        dataset=datasets,
        plot=["scatterplot","heatmap"],
        resolution = config["porec_resolutions"]),
        expand("outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{exp}_{resolution}.insulation.tsv",
        dataset=datasets,
        exp=config["expansions"],
        resolution=config["insulation_resolutions"]),
        expand("outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation_correlation.png",
        dataset=datasets,
        resolution=config["insulation_resolutions"]),
        expand("outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}.cis{file}",
        dataset=datasets,
        exp=config["expansions"],
        file=[".vecs.tsv",".lam.txt",".bw"]),
        expand("outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{resolution}.eigs_correlation.png",
        dataset=datasets,
        resolution=[100000])


rule pairs_to_cooler:
    input:
        fai = f"{config['ref']}.fai",
        pairs = "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{exp}.pairs.gz"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool"
    wildcard_constraints:
        resolution = r"\d+" # = only numbers
    log:
        "logs/pairs_to_cooler.{dataset}_{exp}_{resolution}.log"
    conda:
        "../env/cooler.yml"
    threads: 2
    shell:
        """
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            {input.fai}:{wildcards.resolution} \
            {input.pairs} \
            {output.cool} \
            >{log} 2>&1
        """

# use cooler zoomify to make lower resolutions out of one high res .cool
# how well does this work when the 1000 bp is super sparse? i am using the 5kb instead
# mcool is currently not used anywhere

rule merge_mcools:
    input:
        expand("outputs/pairs_files_T2T/{{dataset}}/cooler/{{dataset}}_{{exp}}_{resolution}.cool", 
        resolution=config.get("min_bin_width", "5000"))
    output:
        mcool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}.mcool"
    log:
        "logs/merge_mcools.{dataset}_{exp}.log"
    params:
        resolutions = config.get("cooler_resolutions", "5000,10000,25000,40000,100000,250000"),
        prefix = lambda wc: wc.dataset
    conda:
        "../env/cooler.yml"
    threads: 2
    shell:
        """
        cooler zoomify \
            -r {params.resolutions} \
            -o {output.mcool} \
            {input}\
            >{log} 2>&1
        """

# use cooler balance instead of hiccorrectmatrix to stay within cooler system
# standard settings
# does ice
# adds weight column to cooler instead of creating new file but i am putting the balanced matrices in a new file

rule cooler_balance:
    input: 
        "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool"
    output: 
        "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_balanced.cool"
    log:
        "logs/cooler_balance.{dataset}_{exp}_{resolution}.log"
    conda:
        "../env/cooler.yml"
    shell:
        """
        cp {input} {output}
        cooler balance \
        --cis-only \
        --ignore-diags 2 \
        --mad-max 5 \
        --min-nnz 10 \
        --tol 1e-05 \
        --max-iters 500 \
        {output} \
        >{log} 2>&1
        """

#A linear correlation between the Pore-C adjacent and non-adjacent contact 
#matrices measured three different similarity metrics:
#(1) raw contact matrices, 
#(2) compartmental eigenvector scores
#(3) insulation scores

# (1) raw contact matrices (5kb, 10kb, 25kb, 40kb)
# hicCorrelate is better used on uncorrected matrices
# here, hicExplorer is used bc no matrix correction has been applied to the raw cooler, 
# if wish to use with corrected matrix, must use hicCorrectmatrix to generate new entries for each bin, not add a weight column

rule correlate_matrices:
    input:
        cool = expand("outputs/pairs_files_T2T/{{dataset}}/cooler/{{dataset}}_{exp}_{{resolution}}.cool",
        exp=config["expansions"])
    output:
        scatterplot = "outputs/pairs_files_T2T/{dataset}/qc/correlation_scatterplot_{resolution}.png",
        heatmap = "outputs/pairs_files_T2T/{dataset}/qc/correlation_heatmap_{resolution}.png"
    params:
        range = "5000:500000" # consider contacts in this range. may too small
    log:
        "logs/correlate_matrices.{dataset}_{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrelate \
        --matrices {input.cool} \
        --method=pearson \
        --log1p \
        --range {params.range} \
        --outFileNameHeatmap {output.heatmap} \
        --outFileNameScatter {output.scatterplot} \
        >{log} 2>&1
        """

# (2) compartmental eigenvector scores (a/b compartments, on 100 kb matrix)
# 1. create phasing track (based on gc content of ref)

rule eigs_phasing_track:
    input: 
        cool = "outputs/pairs_files_T2T/T2T04.hp1/cooler/T2T04.hp1_adj_100000.cool",
        fasta = config['ref']
    output:
        bins = "data/bins.100000.tsv",
        gc = "data/gc.100000.tsv"
    log:
        "logs/eigs_phasing_track.log"
    conda:
        "../env/cooltools.yml"
    shell:
        """
        cooler dump --header -t bins {input.cool} | cut -f1-3 > {output.bins} 2> {log}
        cooltools genome gc {output.bins} {input.fasta} > {output.gc} 2>> {log}
        """

# 2. get eigs of 100 kb matrix, phase based on gc content
# in humans: GC content is useful for phasing because it typically has a strong correlation at the 100kb-1Mb bin level with the eigenvector

rule eigs_cis:
    input:
        gc = "data/gc.100000.tsv",
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_100000_balanced.cool"
    output:
        vecs="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}.cis.vecs.tsv",
        lam="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}.cis.lam.txt",
        bw="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}.cis.bw"
    params:
        prefix = "outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}"
    log:
        "logs/eigs_cis.{dataset}_{exp}.log"
    conda:
        "../env/cooltools.yml"
    shell:
        """
        cooltools eigs-cis \
        --phasing-track {input.gc}::GC \
        --clr-weight-name "weight" \
        --verbose \
        --bigwig \
        --out-prefix {params.prefix} \
        {input.cool} \
        >{log} 2>&1
        """

# correlate eigenvectors (eig1) across two matrices

rule eigs_correlation:
    input:
        adj_eigs= "outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_adj.cis.vecs.tsv",
        nonadj_eigs = "outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_nonadj.cis.vecs.tsv",
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{resolution}.eigs_correlation.png",
    log: 
        "logs/eigs_correlation_{dataset}_{resolution}.log"
    conda: 
        "../env/r_plotting.yml"
    shell:
        """
        Rscript scripts/plot_compartments.R \
        --adj_eigs {input.adj_eigs} \
        --nonadj_eigs {input.nonadj_eigs} \
        --plot {output.plot} \
        >{log} 2>&1
        """

# (3) insulation score
# 1. calculate insulation for 25 and 40 kb
# Li and Otsu automated thresholding criteria borrowed from the image processing field
rule insulation_score:
    input: 
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_balanced.cool"
    output:
        insu = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{exp}_{resolution}.insulation.tsv"
    params:
        window = lambda wildcards: " ".join([str(mult * int(wildcards.resolution)) for mult in config['insu_window_multipliers']])
    log:
        "logs/cooltools_insulation.{dataset}_{exp}_{resolution}.log"
    conda:
        "../env/cooltools.yml"
    shell:
        """
        cooltools insulation \
        --threshold Li \
        --clr-weight-name "weight" \
        --output {output.insu} \
        --bigwig \
        --verbose \
        {input.cool} {params.window} \
        >{log} 2>&1
        """

# 2. plot insulation correlation for a. all bins, b. bins denoted boundaries in nonadj contact map
rule insulation_correlation:
    input:
        adj_insu = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_adj_{resolution}.insulation.tsv",
        nonadj_insu = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_nonadj_{resolution}.insulation.tsv",
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation_correlation.png",
        others = expand("outputs/pairs_files_T2T/{{dataset}}/insulation/{{dataset}}_{{resolution}}.insulation_correlation_{boundary}.png",
        boundary=["3x_boundary_nonadj", "5x_boundary_nonadj", "10x_boundary_nonadj", "25x_boundary_nonadj"])
    params:
        other_plots = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation_correlation_"
    log: 
        "logs/insulation_correlation_{dataset}_{resolution}.log"
    conda: 
        "../env/r_plotting.yml"
    shell:
        """
        Rscript scripts/plot_insulation_corr.R \
        --adj_insu {input.adj_insu} \
        --nonadj_insu {input.nonadj_insu} \
        --prefixes {params.other_plots} \
        --plot {output.plot} \
        >{log} 2>&1
        """

