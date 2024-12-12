from pathlib import Path
import math

def get_bams_from_pod5_split(wildcards):
    wildcards_clean = {'dataset': wildcards.dataset.replace(".subset", "")}
    dataset_split = checkpoints.pod5_split.get(**wildcards_clean).output["dataset_split"]
    bamfiles = [str(x.with_suffix(".bam")) for x in Path(dataset_split).glob("*.pod5")]
    if ".subset" in wildcards.dataset:
        middle = math.floor(len(bamfiles)/2)
        bamfiles = bamfiles[middle:middle+config['process_subset_chunks']]
        bamfiles = [x.replace('.subset', "") for x in bamfiles]
    return bamfiles


def get_path_for_dataset_folder(wc):
    
    # Get basename of folder
    folder_base = wc.dataset
    
    # Search through datasets dictionary for matching folder paths
    matching_paths = []
    for dataset_name, dataset_types in datasets.items():
        for data_type, paths in dataset_types.items():
            # Handle both single path and list of paths
            if isinstance(paths, str):
                paths = [paths]
            
            for path in paths:
                if os.path.isdir(path) and folder_base in path:
                    matching_paths.append(path)
            
            # Print warning if multiple paths are found
            if len(matching_paths) > 1:
                print(f"Warning: Multiple paths found for dataset {dataset_name} and type {data_type}: {matching_paths}")
    #print(matching_paths)                
    return matching_paths

rule dorado_sup:
    input:
        pod5=get_path_for_dataset_folder,
    output:
        done = "data/basecalled/SUP/{dataset,[^.]+(?!\.bam$)}/dorado_sup.done"
    log:
        "logs/dorado_duplex_{dataset}.log",
    resources:
        queue="gpu_srv010,gpu_srv019",
        gpus=2,
    threads: 32
    priority: 3
    params:
        dorado=config["dorado"],
        model=config["model_auto"],
    shell:
        """
        {params.dorado} basecall \
            {params.model} \
            {input.pod5} \
            --models-directory dorado_models \
            --recursive \
            --trim all \
            --output-dir $(dirname {output.done}) \
            2> {log} 2>&1
        touch {output.done}
         """

rule rename_dorado_output:
    input:
        folder="data/basecalled/SUP/{dataset}/dorado_sup.done",
    output:
        bam="data/basecalled/SUP/{dataset}/{dataset}.sup.unmapped.bam",
    log:
        "logs/rename_dorado_output_{dataset}.log",
    shell:
        """
        mv -v $(find {input.folder} -name "calls_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_T[0-9][0-9]-[0-9][0-9]-[0-9][0-9].bam") {output.bam} \
        > {log} 2>&1
        """



rule dorado_summary:
    input:
        folder="data/basecalled/{path}",
    output:
        summary="data/basecalled/{path}/sequencing_summary.txt",
    params:
        dorado=config["dorado"],
    log:
        "logs/dorado_summary_{path}.log",
    shell:
        """
        {params.dorado} summary \
            {input.folder} \
            > {output.summary} \
            2> {log}
        """

"""
Rules using pod5 to preprocess the raw data (extract channel
information, split into channel-specific POD5 files).
"""

rule pod5_view:
    input:
        dataset="{dataset}",
    output:
        summary=temp("{dataset}.summary.tsv"),
    log:
        "logs/{dataset}_pod5_view.log",
    conda:
        "../envs/pod5.yaml"
    priority: 1
    threads: 4
    shell:
        """
        pod5 \
            view \
            --recursive \
            --threads {threads} \
            --include "read_id,channel" \
            --output {output.summary} \
            {input.dataset} \
            > {log} 2>&1
        """

rule process_summary:
    input:
        summary="{dataset}.summary.tsv",
    output:
        summary=temp("{dataset}.summary.proc.tsv"),
    log:
        "logs/{dataset}_process_summary.log",
    params:
        split_number=config["split_number"],
    priority: 1
    run:
        import pandas as pd
        df = pd.read_csv(input["summary"], sep="\t")
        df["channel_mod"] = df["channel"].mod(params["split_number"])
        df.to_csv(output["summary"], sep="\t", index=False)


checkpoint pod5_split:
    input:
        dataset="{dataset}",
        summary="{dataset}.summary.proc.tsv",
    output:
        dataset_split=directory(temp("{dataset}_split_by_channel")),
    log:
        "logs/{dataset}_pod5_split.log",
    conda:
        "../envs/pod5.yaml"
    threads: 20
    priority: 2
    shell:
        """
        pod5 \
            subset \
            {input.dataset} \
            --recursive \
            --threads {threads} \
            --output {output.dataset_split} \
            --summary {input.summary} \
            --columns channel_mod \
            --template "channel-{{channel_mod}}.pod5" \
            > {log} 2>&1
        """



"""
Rules using dorado to basecall POD5 data.
"""

# We still need to manually reduce chunksize and batchsize because of a unsolved bug in dorado
# https://github.com/nanoporetech/dorado/issues/850
# https://github.com/nanoporetech/dorado/issues/842
# https://github.com/nanoporetech/dorado/issues/1098

rule dorado_duplex:
    input:
        pod5="{dataset}_split_by_channel/channel-{channel}.pod5",
    output:
        bam=temp("{dataset}_split_by_channel/channel-{channel}.bam"),
    log:
        "logs/dorado_duplex_{dataset}_{channel}.log",
    resources:
        queue="gpu_srv010,gpu_srv019",
        gpus=1,
    threads: 32
    priority: 3
    params:
        dorado=config["dorado"],
        model=config["model_auto_duplex"],
    run:
        with get_gpu_id() as gid:  # Check for unused GPU
            params.cuda_device = f"cuda:{gid}"
            shell(
                "CUDA_LAUNCH_BLOCKING=1 \
                {params.dorado} duplex \
                {params.model} \
                {input.pod5} \
                --models-directory dorado_models \
                --device '{params.cuda_device}' \
                --chunksize 9996 \
                --batchsize 500 \
                > {output.bam} \
                2> {log} \
            "
            )


rule concat_bam_filelist:
    input:
        bams=get_bams_from_pod5_split,
    output:
        dataset_split_bams=temp("{dataset}.bam.fofn"),
    run:
        with open(output["dataset_split_bams"], "w") as out:
            for bam in input["bams"]:
                out.write(bam + "\n")


rule concat_bam:
    input:
        dataset_split_bams="{dataset}.bam.fofn",
        bams=get_bams_from_pod5_split,
    output:
        bam="{dataset}.duplex.unmapped.bam",
    log:
        "logs/concat_bam_{dataset}.log",
    conda:
        "../envs/samtools.yaml"
    priority: 4
    threads: 4
    shell:
        """
        samtools \
            cat \
            --threads {threads} \
            -o {output.bam} \
            -b {input.dataset_split_bams} \
            > {log} 2>&1
        """

rule report_duplex_statistics:
    input:
        bam="{dataset}.duplex.unmapped.bam",
    output:
        statistics="{dataset}.duplex.stat.txt",
    log:
        "logs/report_duplex_statistics_{dataset}.log",
    conda:
        "../env/samtools.yml"
    threads: 1
    params:
        script=workflow.source_path("../scripts/bam-duplex-rate.awk"),
    priority: 5
    shell:
        """
        samtools view {input.bam} 2> {log} | \
            awk -f {params.script} >> {output.statistics} 2>> {log}
        """

rule extract_duplex:
    input:
        bam="data/basecalled/Duplex/{dataset}/{dataset}.duplex.unmapped.bam",
        stats = "data/basecalled/Duplex/{dataset}/{dataset}.duplex.stat.txt",
    output:
        bam="data/basecalled/Duplex/{dataset}/{dataset}.duplexonly.unmapped.bam"
    log:
        "logs/extract_duplex_{dataset}.log"
    shell:
        """
        samtools view -b -h -d dx:1 {input.bam} > {output.bam} 2> {log}
        """
