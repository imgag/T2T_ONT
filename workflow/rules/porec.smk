rule all_porec:
    input:
        "analysis_other/porec/T2T04/pairs/T2T04.pairs.gz",
        "analysis_other/porec/T2T04/pairs/T2T04.pairs.stats.html",

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