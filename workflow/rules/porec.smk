rule all_porec:
    input:
        "analysis_other/porec/T2T04/pairs/T2T04.pairs.gz",
        "analysis_other/porec/T2T04/pairs/T2T04.pairs.stats.html",
        "analysis_other/porec/T2T04/cooler/T2T04.mcool",
        expand("analysis_other/porec/T2T04/cooler/T2T04_{res}_corrected.cool", res = config['porec_resolutions']),
        "analysis_other/porec/T2T04/qc/T2T04_10000_diagnostic.png",
        "analysis_other/porec/T2T04/qc/plot_vs_counts_10000.png",
        "analysis_other/porec/T2T04/qc/matrix_correlation_heatmap_10000.png"


def get_all_porec_runs(wc):
    folders = datasets[wc.dataset].get("POREC", [])
    # If folders is not a list, make it one
    if not isinstance(folders, list):
        folders = [folders] if folders else []
    
    run_ids = []
    for folder_path in folders:
        # Extract the basename (run ID) from the folder path
        run_id = os.path.basename(folder_path)
        # Look up the run_id in unique_datasets["porec"] to verify it exists
        if run_id in unique_datasets["porec"]:
            run_ids.append(run_id)

    return run_ids

rule merge_pairs:
    input:
        pairs = lambda wc: expand("analysis_other/wf-pore-c/{run}/pairs/{run}.pairs.gz", run = get_all_porec_runs(wc))
    output:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    log:
        "logs/porec/merge_pairs.{dataset}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools merge \
            --output {output.pairs} \
            {input.pairs} \
            >{log} 2>&1
        """

rule merge_pairs_stats:
    input:
        stats = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    output:
        stats = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.txt"
    log:
        "logs/porec/merge_pairs_stats.{dataset}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools stats \
            --output {output.stats} \
            {input.stats} \
            >{log} 2>&1
        """

rule pairs_stats_report:
    input:
        "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.txt"
    output:
        "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.html"
    log:
        "logs/porec/merge_pairs_stats_report.{dataset}.log"
    params:
        report_script = "workflow/scripts/create_pairs_report.py"
    conda:
        "../env/py_report.yml"
    shell:
        """
        python {params.report_script} \
            {input} {output} \
            >{log} 2>&1
        """


rule pairs_to_cooler:
    input:
        fai = f"{config['ref']}.fai",
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.cool"
    wildcard_constraints:
        resolution = "\d+"
    log:
        "logs/porec/pairs_to_cooler.{dataset}_{resolution}.log"
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

rule merge_mcools:
    input:
        "analysis_other/porec/{dataset}/cooler/{dataset}_1000.cool", 
    output:
        mcool = "analysis_other/porec/{dataset}/cooler/{dataset}.mcool"
    log:
        "logs/porec/merge_mcools.{dataset}.log"
    params:
        resolutions = config.get("cooler_resolutions", "1000,5000,10000,50000,100000"),
        prefix = lambda wc: wc.dataset
    conda:
        "../env/cooler.yml"
    threads: 2
    shell:
        """
        cooler zoomify \
            -r {params.resolutions} \
            -o {output.mcool} \
            {params.prefix}.cool \
            >{log} 2>&1
        """

rule hic_normalize:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.cool"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool"
    log:
        "logs/porec/hic_normalize.{dataset}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicNormalize \
            --matrices {input.cool} \
            --normalize smallest \
            -o {output.cool} \
            >{log} 2>&1
        """

rule hic_diagnostic_plot:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool"
    output:
        plot = "analysis_other/porec/{dataset}/qc/{dataset}_{resolution}_diagnostic.png"
    log:
        "logs/porec/hic_diagnostic.{dataset}.{resolution}.log"
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
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    log:
        "logs/porec/hic_correct.{dataset}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),
        filter_threshold = config.get("filter_threshold", "-2.5 5")
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

rule hic_plot_dist_vs_counts:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        plot = "analysis_other/porec/{dataset}/qc/plot_vs_counts_{resolution}.png"
    log:
        "logs/porec/plot_dist_vs_counts.{dataset}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotDistVsCounts \
            --matrices {input.cool} \
            -o {output.plot} \
            >{log} 2>&1
        """

rule hic_correlate:
    input:
        cool_files = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        heatmap = "analysis_other/porec/{dataset}/qc/matrix_correlation_heatmap_{resolution}.png",
        scatterplot = "analysis_other/porec/{dataset}/qc/matrix_correlation_scatter_{resolution}.png"
    log:
        "logs/porec/hic_correlate.{dataset}.{resolution}.log"
    params:
        range_param = config.get("correlation_range", "20000:500000")
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrelate \
            --log1p \
            --matrices {input.cool_files} \
            --range {params.range_param} \
            -oh {output.heatmap} \
            -os {output.scatterplot} \
            >{log} 2>&1
        """
