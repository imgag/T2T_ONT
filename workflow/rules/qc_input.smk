rule mosdepth:
    input:
        bam="data/mapped/{path}.bam",
    output:
        depth="data/bamstats/{path}/cov.mosdepth.summary.txt",
        dist="data/bamstats/{path}/cov.mosdepth.global.dist.txt",
    threads: 4
    conda:
        "../env/mosdepth.yml"
    log:
        "logs/modepth_{path}.log",
    shell:
        """
        mosdepth \
            -n --fast-mode --by 500 \
            -t {threads}\
            $(dirname {output.depth})/cov \
            {input.bam} \
            >{log} 2>&1
        """


rule bam_qc:
    input:
        bam="data/mapped/{path}.bam",
        ref=config["ref"],
    output:
        tsv="data/bamstats/{path}/bamqc.tsv.gz",
        json="data/bamstats/{path}/bamqc.json.gz",
    log:
        "logs/bam_qc_{path}.log",
    threads: 1
    params:
        alfred="bin/alfred",
    shell:
        """
        {params.alfred} qc \
            -r {input.ref} \
            -o {output.tsv} \
            -j {output.json} \
            -i {input.bam} \
            >{log} 2>&1
        """


rule bamstats:
    input:
        bam="data/mapped/{path}.bam",
    output:
        hist_l="data/bamstats/mapped/{path}/length.hist",
        hist_q="data/bamstats/mapped/{path}/quality.hist",
        hist_a="data/bamstats/mapped/{path}/accuracy.hist",
        hist_c="data/bamstats/mapped/{path}/coverage.hist",
        bamstats="data/bamstats/mapped/{path}/bamstats.txt",
        flagstats="data/bamstats/mapped/{path}/flagstats.txt",
        basecallers="data/bamstats/mapped/{path}/basecallers.txt",
    conda:
        "../env/fastcat.yml"
    log:
        "logs/bamstats.{path}.log",
    threads: 4
    shell:
        """
        tmp=$(mktemp -u)
        bamstats \
            --histograms=$tmp \
            --flagstats=$tmp/$(basename {output.flagstats}) \
            --basecallers=$tmp/$(basename {output.basecallers}) \
            --sample={wildcards.path}\
            --threads {threads} \
            {input.bam} >{output.bamstats} 2>{log}
        cp -r $tmp/* $(dirname {output.bamstats})
        rm -rf $tmp
        """


rule fastcat_qc_only:
    input:
        fq="assembly/input/{dataset}/{dataset}.{type}.fastq.gz",
    output:
        hist_l="assembly/input_qc/{dataset}/{dataset}.{type}/length.hist",
        hist_q="assembly/input_qc/{dataset}/{dataset}.{type}/quality.hist",
        stat_reads="assembly/input_qc/{dataset}/{dataset}.{type}/read_stats.txt",
    conda:
        "../env/fastcat.yml"
    log:
        "logs/fastcat_qc_only.{dataset}.{type}.log",
    shell:
        """
        tmp=$(mktemp -u)
        fastcat \
            --histograms=$tmp\
            --read=$tmp/$(basename {output.stat_reads}) \
            --sample={wildcards.dataset}_{wildcards.type} \
            {input} >/dev/null 2>{log}
        cp -r $tmp/* $(dirname {output.stat_reads})
        rm -rf $tmp
        """


rule process_dorado_summary:
    input:
        "data/basecalled/{type}/{dataset}/{dataset}.{model}.sequencing_summary.txt",
    output:
        "data/bamstats/basecalled/{type}/{dataset}/{dataset}.{model}.sequencing_summary.processed.txt",
    conda:
        "../env/R.yml"
    log:
        "logs/process_dorado_summary.{dataset}_{model}_{type}.txt"
    shell:
        """
        Rscript workflow/scripts/07_dorado_summary_stats.R {input} {output} > {log} 2>&1
        """
