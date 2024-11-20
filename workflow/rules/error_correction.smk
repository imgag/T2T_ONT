rule dorado_trim:
    input:
        fastq=lambda wc: datasets[wc.dataset].get("HQ_herro", ""),
    output:
        fastq="data/corrected/{dataset}/{file}.trimmed.fastq",
    params:
        dorado=config["dorado"],
    threads: 60
    benchmark:
        "runtimes/dorado_trim_{dataset}_{file}.txt"
    log:
        "logs/dorado_trim_{dataset}_{file}.log",
    shell:
        """
        {params.dorado} trim \
            --threads {threads} \
            --emit-fastq \
            {input.fastq} > {output.fastq} \
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
    threads: 60
    benchmark:
        "runtimes/dorado_correct_mapping_{dataset}_{file}.txt"
    log:
        "logs/dorado_correct_mapping_{dataset}_{file}.log",
    shell:
        """
        {params.dorado} correct \
            --to-paf \
            --threads {threads} \
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
    threads: 6
    benchmark:
        "runtimes/dorado_correct_inference_{dataset}_{file}.txt"
    log:
        "logs/dorado_correct_inference_{dataset}_{file}.log",
    resources:
        gpu=2,
    shell:
        """
        {params.dorado} correct \
            --from-paf {input.paf} \
            --threads {threads} \
            --model-path {params.herro_model} \
            --device 'cuda:all' \
            {input.fastq} > {output.fa} \
            2> {log}
       """
