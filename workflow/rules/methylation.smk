rule all_methylation:
    input:
        expand("analysis_other/mod/{sample}/{sample}.modkit_summary.tsv", sample = "T2T01"),
        expand("analysis_other/mod/{sample}/{sample}.modifications.bed", sample = "T2T01")


rule modkit_summary:
    input:
        "data/mapped/{sample}/{sample}.UL.asm.bam"
    output:
        "analysis_other/mod/{sample}/{sample}.modkit_summary.tsv"
    log:
        "logs/methylation/{sample}.modkit_summary.log"
    params:
        modkit = config['modkit'],
        mod_num_reads = config['mod_num_reads']
    threads:
        8
    shell:
        """
        {params.modkit} summary \
            --threads {threads} \
            --log-filepath {log} \
            --num-reads {params.mod_num_reads} \
            {input} > {output} \
            2>>{log}
        """


rule modkit_pileup:
    input:
        "data/mapped/{sample}/{sample}.UL.asm.bam"
    output:
        "analysis_other/mod/{sample}/{sample}.modifications.bed"
    log:
        "logs/methylation/{sample}.modkit_summary.log"
    params:
        modkit = config['modkit'],
        mod_num_reads = config['mod_num_reads']
    threads:
        8
    shell:
        """
        {params.modkit} pileup \
            --threads {threads} \
            --num-reads {params.mod_num_reads} \
            --log-filepath {log} \
            {input} {output} \
            2>>{log}
        """