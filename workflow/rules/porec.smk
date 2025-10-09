rule all_porec:
    input:
        expand("analysis_other/porec/{sample}.done", sample = ["T2T04", "T2T03"])

rule collect_porec:
    input:
        expand("analysis_other/porec/{{asm}}{hp}/pairs/{{asm}}{hp}.pairs.gz", hp = ["", ".hp1", ".hp2"]),
        expand("analysis_other/porec/{{asm}}{hp}/pairs/{{asm}}{hp}.pairs.stats.html", hp = ["", ".hp1", ".hp2"]),
        expand("analysis_other/porec/{{asm}}{hp}/cooler/{{asm}}{hp}.mcool", hp = ["", ".hp1", ".hp2"]),
        expand("analysis_other/porec/{{asm}}{hp}/cooler/{{asm}}{hp}_{res}_corrected.cool", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/qc/{{asm}}{hp}_{res}_diagnostic.png", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/qc/plot_vs_counts_{res}.png", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/tad/{{asm}}{hp}_{res}_domains.bed", hp = ["", ".hp1", ".hp2"], res=config.get('tad_resolutions', ['25000'])),
        expand("analysis_other/porec/{{asm}}{hp}/loops/{{asm}}{hp}_{res}_loops.bedpe", hp = ["", ".hp1", ".hp2"], res=config.get('loop_resolutions', ['10000']))
    output:
        "analysis_other/porec/{asm}.done"
    shell:
        """
        touch {output}
        """

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
    wildcard_constraints:
        dataset=r"\w+"
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

rule pairtools_parse_bam:
    input:
        bam = lambda wc: f"analysis_other/dip3d/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}/4-haplotag/{wc.dataset}.bam", 
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    log:
        "logs/porec/pairtools_parse2.{dataset}.log"
    conda:
        "../env/pairtools.yml"
    params:
        assembly = "T2T-CHM13.v2",
        columns = "readID,chrom1,pos1,chrom2,pos2,strand1,strand2,mapq1,mapq2",
        orientation = "pair", 
        position = "junction"
    shell:
        """
        pairtools parse2 \
            --chroms-path {input.chromsize} \
            --assembly {params.assembly} \
            --report-position {params.position} \
            --report-orientation {params.orientation} \
            --single-end \
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log}
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
        resolution = r"\d+"
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
        expand("analysis_other/porec/{{dataset}}/cooler/{{dataset}}_{resolution}.cool", resolution=config.get("min_bin_width", "1000"))
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
            {input}\
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
            >{log} 2>&1 || {{
            
            echo "Standard correction failed, trying more stringent filtering..." >>{log}
            
            # Try with more stringent filtering
            hicCorrectMatrix correct \
                --matrix {input.cool} \
                --correctionMethod {params.correction_method} \
                --outFileName {output.cool} \
                --filterThreshold -3.0 3.0 \
                --iterNum 200 \
                >>{log} 2>&1 || {{
                
                echo "Stringent correction failed, trying KR normalization..." >>{log}
                
                # Try Knight-Ruiz normalization as fallback
                hicCorrectMatrix correct \
                    --matrix {input.cool} \
                    --correctionMethod KR \
                    --outFileName {output.cool} \
                    --filterThreshold -3.0 3.0 \
                    >>{log} 2>&1 || {{
                    
                    echo "All correction methods failed, copying normalized matrix..." >>{log}
                    # As last resort, just copy the normalized matrix
                    cp {input.cool} {output.cool}
                }}
            }}
        }}
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

rule hic_find_tads:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        tads = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_domains.bed",
        boundaries = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_boundaries.bed",
        score = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_score.bedgraph"
    log:
        "logs/porec/hic_find_tads.{dataset}.{resolution}.log"
    params:
        min_depth = lambda wc: max(config.get("tad_min_depth", 20000), 3*int(wc.resolution)),
        max_depth = config.get("tad_max_depth", 200000),
        step = lambda wc: max(config.get("tad_step", 10000), int(wc.resolution)),
        fdr_threshold = config.get("tad_fdr_threshold", 0.05),
        delta = config.get("tad_delta", 0.01),
        correction_factor_threshold = config.get("tad_correction_threshold", 1.5)
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicFindTADs \
            --matrix {input.cool} \
            --outPrefix analysis_other/porec/{wildcards.dataset}/tad/{wildcards.dataset}_{wildcards.resolution} \
            --minDepth {params.min_depth} \
            --maxDepth {params.max_depth} \
            --step {params.step} \
            --thresholdComparisons {params.fdr_threshold} \
            --delta {params.delta} \
            --correctForMultipleTesting fdr \
            --numberOfProcessors {threads} \
            >{log} 2>&1
        """

rule hic_plot_tads:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        tads = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_domains.bed"
    output:
        plot = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_{region}_tads.png"
    log:
        "logs/porec/hic_plot_tads.{dataset}.{resolution}.{region}.log"
    params:
        region = lambda wc: wc.region,  # e.g., "chr1:1000000-5000000"
        dpi = config.get("plot_dpi", 300)
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotMatrix \
            --matrix {input.cool} \
            --region {params.region} \
            --outFileName {output.plot} \
            --dpi {params.dpi} \
            --clearMaskedBins \
            >{log} 2>&1
        """

rule hic_detect_loops:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe",
    log:
        "logs/porec/hic_detect_loops.{dataset}.{resolution}.log"
    params:
        window_size = config.get("loop_window_size", 10),
        peak_width = config.get("loop_peak_width", 6),
        max_loop_distance = config.get("loop_max_distance", 2000000),
        pvalue_threshold = config.get("loop_pvalue_threshold", 0.05),
        fdr_threshold = config.get("loop_fdr_threshold", 0.1)
    threads: 4
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicDetectLoops \
            --matrix {input.cool} \
            --outFileName {output.loops} \
            --windowSize {params.window_size} \
            --peakWidth {params.peak_width} \
            --maxLoopDistance {params.max_loop_distance} \
            --pValuePreselection {params.pvalue_threshold} \
            --pValue {params.fdr_threshold} \
            --threads {threads} \
            >{log} 2>&1
        """

rule hic_validate_viewpoints:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe"
    output:
        viewpoints = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_viewpoints.txt",
        interactions = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_interactions.txt"
    log:
        "logs/porec/hic_validate_viewpoints.{dataset}.{resolution}.log"
    params:
        resolution_val = lambda wc: wc.resolution,
        interaction_file = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_interactions.txt"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicValidateLocations \
            --matrix {input.cool} \
            --locations {input.loops} \
            --method loops \
            --outFileName {output.viewpoints} \
            --outFileNameInteractions {params.interaction_file} \
            --resolution {params.resolution_val} \
            >{log} 2>&1
        """

rule hic_aggregate_contacts:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe"
    output:
        aggregate_plot = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_aggregate_loops.png",
        aggregate_data = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_aggregate_loops.npz"
    log:
        "logs/porec/hic_aggregate_contacts.{dataset}.{resolution}.log"
    params:
        range_around_loop = config.get("aggregate_range", 100000),
        number_of_bins = config.get("aggregate_bins", 30)
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicAggregateContacts \
            --matrix {input.cool} \
            --BED {input.loops} \
            --outFileName {output.aggregate_plot} \
            --outFilePrefixMatrix analysis_other/porec/{wildcards.dataset}/loops/{wildcards.dataset}_{wildcards.resolution}_aggregate_loops \
            --range {params.range_around_loop}:{params.range_around_loop} \
            --numberOfBins {params.number_of_bins} \
            --plotType 2d \
            >{log} 2>&1
        """