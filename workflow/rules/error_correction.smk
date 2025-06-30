def get_dorado_trim_input(wc):
    """Retuns fastq if from published dataset, else SUP basecalled"""
    if [wc.dataset] == "published":
        hts_file = datasets[wc.dataset].get("HQ_herro", "")
    else:
        hts_file = f"data/basecalled/sup/{wc.file}/{wc.file}.sup.unmapped.bam"
    return hts_file

# We need fastq as input for dorado correct. Might as well trim it during conversion.
rule dorado_trim:
    input:
        get_dorado_trim_input,
    output:
        fastq="data/corrected/{dataset}/{file}.trimmed.fastq",
    params:
        dorado=config["dorado"],
    threads: 5
    benchmark:
        "runtimes/dorado_trim/{dataset}_{file}.txt"
    log:
        "logs/dorado_trim/{dataset}_{file}.log",
    shell:
        """
        {params.dorado} trim \
            --threads 20 \
            --emit-fastq \
            {input} > {output.fastq} \
            2> {log}
        """


rule dorado_correct_mapping:
    input:
        fastq=rules.dorado_trim.output.fastq,
    output:
        paf="data/corrected/{dataset}/{file}.ovl.paf",
    params:
        dorado=config["dorado"],
        herro_model=config["herro_model"],
#    threads: 30 get more done
    threads: 5
    benchmark:
        "runtimes/dorado_correct_mapping/{dataset}_{file}.txt"
    log:
        "logs/dorado_correct_mapping/{dataset}_{file}.log",
    shell:
        """
        {params.dorado} correct \
            --to-paf \
            --threads 30 \
            --model-path {params.herro_model} \
            --device 'cpu' \
            {input.fastq} > {output.paf} \
            2> {log}
        """


rule dorado_correct_inference:
    input:
        fastq=rules.dorado_trim.output.fastq,
        paf=rules.dorado_correct_mapping.output.paf,
    output:
        fa="data/corrected/{dataset}/{file}.corrected.fasta",
    params:
        dorado=config["dorado"],
        herro_model=config["herro_model"],
    threads: 1
    benchmark:
        "runtimes/dorado_correct_inference/{dataset}_{file}.txt"
    log:
        "logs/dorado_correct_inference/{dataset}_{file}.log",
    resources:
        queue="gpu_srv010", 
        gpu=1,
    run:
        with get_gpu_id() as gid:  # Check for unused GPU
            params.cuda_device = f"cuda:{gid}"
            shell(
                "{params.dorado} correct \
                --from-paf {input.paf} \
                --threads {threads} \
                --model-path {params.herro_model} \
                --device '{params.cuda_device}' \
                --index-size 4G \
                {input.fastq} > {output.fa} \
                2> {log} \
            "
            )
