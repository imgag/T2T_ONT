#ice correction with thresholds taken from diagnostic plot log output

rule hic_correct_matrix:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool"
        diagnostic = "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_corrected.cool"
    log:
        "logs/hic_correct.{dataset}_{exp}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),
        #filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        madscore=$(grep "mad threshold" {input.diagnostic} | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g')
        upper=$(echo -3*$madscore | bc)
        thresholds=$(echo $madscore " " $upper)

        hicCorrectMatrix correct \
        --filterThreshold $thresholds 
        --matrix {input.cool} \
        --correctionMethod {params.correction_method} \
        --outFileName {output.cool} \
        >{log} 2>&1
        """
# old:
# filter_threshold need to be set individually. but not quite feasible - what is a very conservative choice?
rule hic_correct_matrix:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool",
        diagnostic = "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_corrected.cool"
    log:
        "logs/hic_correct.{dataset}_{exp}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),
        filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrectMatrix correct \
            --matrix {input.cool} \
            --correctionMethod {params.correction_method} \
            --outFileName {output.cool} \
            --filterThreshold {params.filter_threshold} \
            >{log} 2>&1
        """
        
rule expected_cis:
    shell:
        """
        cooltools expected-cis \
        --outout {output.expected} \
        --smooth --aggregate-smoothed \
        {input.cool} \
        >{log} 2>&1
        """

# 3. Correlate all matrices?
rule correlate_matrices:
    input:
        cool = expand("outputs/pairs_files_T2T/{dataset}{hp}/cooler/{dataset}{hp}_{resolution}.cool",
        dataset = finished_samples,
        hp = ["", ".hp1", ".hp2"],
        resolution = 100000)
    output:
        scatterplot = "outputs/pairs_files_T2T/combined_plots/correlation_scatterplot_{resolution}.png",
        heatmap = "outputs/pairs_files_T2T/combined_plots/correlation_heatmap_{resolution}.png"
    params:
        range = "5000:5000000" # consider contacts in this range. how to decide?
    log:
        "logs/porec/correlate_matrices.{dataset}_{resolution}.log"
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


#cooltools dots 
#--nproc 6 -o 'outputs/test.dots.10000.tsv' 
#--view "data/view_hg38.tsv" 
#$cool_file::resolutions/10000 
#outputs/test.expected.cis.10000.tsv

datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule fanc_triangle_plot:
    input:
        hic=
        insulation=
    output:
        plot=
    params:
        resolution=
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.{dataset}_{exp}_{resolution}.log"
    shell:
        """
        fancplot 
        -o {output.plot} chr18:18mb-28mb \
        -p {input.hic}@{params.resolution} \
        -m 4000000 \
        -vmin 0 -vmax 50 \
        -p line {input.insulation}
        """

-------------------------------------
rule fanc_plot_whole_chr_eigs:
    input:
        hic= "outputs/hic_files_T2T/{dataset}/hic/{dataset}_nonadj.hic",
        eigs="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_nonadj.cis.bw"
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/whole_chr/{i}_{{dataset}}_{chrom}_None-None.pdf",
        zip,
        i=indices_chrom,
        chrom = chroms)
    params:
        chroms=config['chroms_fanc'],
        resolution=250000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/whole_chr",
        vmax=200
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.whole_chr.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.chroms} \
        -p square \
        -vmin 0 -vmax {params.vmax} \
        {input.hic}@{params.resolution}@KR \
        -p line -f -c black --fix-chromosome {input.eigs} \
        >{log} 2>&1
        """


rule fanc_plot_imprinted_genes:
    input:
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic", 
        insulation1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj_25000.insulation.tsv.125000.bw", 
        insulation2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj_25000.insulation.tsv.125000.bw",
        
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/imprinted_genes/{i}_{{dataset}}_{region}.pdf",
        zip,
        i=indices_genes,
        region = imprinted_genes_filename)
    params:
        regions=config['regions_fanc'],
        titles = config['regions_titles_fanc'],
        resolution=10000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes",
        vmax=35
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.imprinted_genes.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.regions} \
        -p mirror \
        -uvmin 0 -uvmax {params.vmax} \
        -lvmin 0 -lvmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -l hp1 hp2 -c black red --fix-chromosome {input.insulation1} {input.insulation2} \
        >{log} 2>&1
        """

# -c for color, -l for log, -d max distance (pyramid)

rule fanc_plot_insulation:
    input:
        hic= "outputs/hic_files_T2T/{dataset}/hic/{dataset}_nonadj.hic",
        insulation="outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_nonadj_25000.insulation.tsv.125000.bw",
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/insulation_chr_arms/{i}_{{dataset}}_{chrom}.pdf",
        zip,
        i=indices_chrom_arms,
        chrom = chrom_arms)
    params:
        chroms=config['chrom_arms_fanc'],
        resolution=25000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/insulation_chr_arms",
        vmax=100
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.chrom_arms_insu.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.chroms} \
        -p triangular \
        -vmin 0 -vmax {params.vmax} \
        {input.hic}@{params.resolution}@KR \
        -p line --fix-chromosome {input.insulation} \
        >{log} 2>&1
        """

   #-p gene --fix-chromosome {input.genes} \
