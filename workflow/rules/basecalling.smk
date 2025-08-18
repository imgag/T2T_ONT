from pathlib import Path
import math


def get_bams_from_pod5_split(wildcards):
    wildcards_clean = {"dataset": wildcards.dataset.replace(".subset", "")}
    dataset_split = checkpoints.pod5_split.get(**wildcards_clean).output[
        "dataset_split"
    ]
    bamfiles = [str(x.with_suffix(".bam")) for x in Path(dataset_split).glob("*.pod5")]
    if ".subset" in wildcards.dataset:
        middle = math.floor(len(bamfiles) / 2)
        bamfiles = bamfiles[middle : middle + config["process_subset_chunks"]]
        bamfiles = [x.replace(".subset", "") for x in bamfiles]
    return bamfiles


rule dorado:
    input:
        pod5=lambda wc: unique_datasets[wc.type][os.path.basename(wc.dataset)],
    output:
        done=r"data/basecalled/{type}/{dataset,[^.]+(?!\.bam$)}/dorado_{type}.done",
    log:
        "logs/dorado/{type}_{dataset}.log"
    resources:
        queue="gpu_srv010", 
        gpu=2
    benchmark:
        "runtimes/dorado/{type}/{dataset}.{type}.txt"
    threads: 2
    priority: 3
    params:
        dorado=config["dorado"],
        model=lambda wc: config["dorado_model"][wc.type],
        models_directory=config["models_directory"],
        model_mod = lambda wc: f"--modified-bases-models {config["dorado_model"]["mod"]}" if wc.type == "sup" else ""
    shell:
        """
        {params.dorado} basecaller \
            {params.model} \
            {input.pod5} \
            --models-directory {params.models_directory} {params.model_mod} \
            --recursive \
            --trim all \
            --output-dir $(dirname {output.done}) \
            2> {log}
        touch {output.done}
        """

rule rename_dorado_output:
    input:
        folder="data/basecalled/{type}/{dataset}/dorado_{type}.done",
    output:
        bam="data/basecalled/{type}/{dataset}/{dataset}.{type}.unmapped.bam",
    localrule:
        True
    log:
        "logs/rename_dorado_output/{dataset}_{type}.log",
    shell:
        """
        mv -v $(find $(dirname {input.folder}) -name "calls_[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]_T[0-9][0-9]-[0-9][0-9]-[0-9][0-9].bam") {output.bam} \
        > {log} 2>&1
        """


rule dorado_summary:
    input:
        folder="data/basecalled/{type}/{dataset}/{dataset}.{model}.unmapped.bam",
    output:
        summary="data/basecalled/{type}/{dataset}/{dataset}.{model}.sequencing_summary.txt",
    localrule:
        True
    params:
        dorado=config["dorado"],
    benchmark:
        "runtimes/dorado_summary/{dataset}.dorado_summary.{model}.{type}.txt"
    log:
        "logs/dorado_summary/{dataset}_{model}_{type}.log",
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
        pod5=lambda wc: unique_datasets["duplex"][os.path.basename(wc.dataset)],
    output:
        summary=temp(r"{dataset,[^.]+(?!\.bam$)}.summary.tsv"),
    log:
        "logs/pod5_view/{dataset}.log",
    conda:
        "../env/pod5.yml"
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
            {input.pod5} \
            > {log} 2>&1
        """


rule process_summary:
    input:
        summary="{dataset}.summary.tsv",
    output:
        summary=temp("{dataset}.summary.proc.tsv"),
    localrule:
        True
    log:
        "logs/process_summary/{dataset}.log"
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
        dataset=lambda wc: unique_datasets["duplex"][os.path.basename(wc.dataset)],
        summary="{dataset}.summary.proc.tsv",
    output:
        dataset_split=directory(temp("{dataset}_split_by_channel")),
    log:
        "logs/pod5_split/{dataset}.log"
    conda:
        "../env/pod5.yml"
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
        "logs/dorado_duplex/{dataset}_{channel}.log"
    resources:
        queue=config["gpu_queues"],
        gpu=1
    threads: 1
    priority: 3
    params:
        dorado=config["dorado"],
        model=lambda wc: config["dorado_model"]["duplex"],
        models_directory=config["models_directory"],
    run:
        with get_gpu_id() as gid:  # Check for unused GPU
            params.cuda_device = f"cuda:{gid}"
            shell(
                r"CUDA_LAUNCH_BLOCKING=1 \
                {params.dorado} duplex \
                    {params.model} \
                    {input.pod5} \
                    --models-directory {params.models_directory} \
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
    localrule:
        True
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
        "logs/concat_bam/{dataset}.log",
    conda:
        "../env/samtools.yml"
    priority: 4
    threads: 4
    shell:
        """
        samtools \
            cat \
            -o {output.bam} \
            -b {input.dataset_split_bams} \
            > {log} 2>&1
        """


rule report_duplex_statistics:
    input:
        bam="{dataset}.duplex.unmapped.bam",
    output:
        statistics="{dataset}.duplex.stat.txt",
    localrule:
        True
    log:
        "logs/report_duplex_statistics/{dataset}.log",
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
        stats="data/basecalled/Duplex/{dataset}/{dataset}.duplex.stat.txt",
    localrule:
        True
    output:
        bam="data/basecalled/Duplex/{dataset}/{dataset}.duplexonly.unmapped.bam",
    conda:
        "../env/samtools.yml"
    log:
        "logs/extract_duplex/{dataset}.log",
    shell:
        """
        samtools view -b -h -d dx:1 {input.bam} > {output.bam} 2> {log}
        """
