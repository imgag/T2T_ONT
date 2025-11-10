datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_hicExplorer:
    input:
       #expand("outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}.mcool",
       #dataset=datasets,
       #exp=config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/qc/{dataset}_{exp}_{resolution}_diagnostic.png",
        dataset=datasets,
        resolution = config["porec_resolutions"],
        exp=config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/qc/plot_vs_counts_{resolution}.png",
        dataset=datasets,
        resolution = config["porec_resolutions"])


rule pairs_to_cooler:
    input:
        fai = config['chromsize_porec'],
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
# how well doe this work when the 1000 bp is super sparse? i am using th 5kb instead

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

# matrix correction 
rule hic_diagnostic_plot:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool"
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/qc/{dataset}_{exp}_{resolution}_diagnostic.png"
    log:
        "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrectMatrix diagnostic_plot \
            --matrix {input.cool} \
            -o {output.plot} \
            >{log} 2>&1
        """

rule hic_correct_matrix:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool",
        diagnostic = "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_corrected.cool",
        lower_file=temp("outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.lower"),
        upper_file=temp("outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.upper")
    log:
        "logs/hic_correct.{dataset}_{exp}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),

        #filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        echo "Start log ....." > {log}
        grep "mad threshold" {input.diagnostic} 2>> {log} | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g' 2>> {log} > {output.lower_file}
        echo -3*$(cat {output.lower_file}) | bc 2>>{log} > {output.upper_file}

        hicCorrectMatrix correct \
        --filterThreshold $(cat {output.lower_file}) $(cat {output.upper_file}) \
        --matrix {input.cool} \
        --correctionMethod {params.correction_method} \
        --outFileName {output.cool} \
        >>{log} 2>&1
        
        """

# compare the matrices

rule hic_plot_dist_vs_counts:
    input:
        cool = expand("outputs/pairs_files_T2T/{{dataset}}/cooler/{{dataset}}_{exp}_{{resolution}}_corrected.cool",
        exp=config["expansions"])
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/qc/plot_vs_counts_{resolution}.png"
    log:
        "logs/plot_dist_vs_counts.{dataset}_{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotDistVsCounts \
            --matrices {input.cool} \
            -o {output.plot} \
            >{log} 2>&1
        """


