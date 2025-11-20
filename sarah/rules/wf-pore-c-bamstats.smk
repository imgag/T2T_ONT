import yaml

with open("data/datasets.yml") as f:
    datasets = yaml.safe_load(f)

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


rule fragments:
    input:
        pairs = lambda wc: expand("../analysis_other/wf-pore-c/{run}/bams/{run}.ns.bam", run = get_all_porec_runs(wc))
    output:
        frag_count = "outputs/"

    
