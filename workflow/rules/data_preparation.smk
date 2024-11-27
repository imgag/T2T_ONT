def update_herro_paths(f, dataset):
    f = os.path.basename(f)
    match = re.search(
        r"(\.fastq|\.fastq.gz|\.bam|\.fasta|\.fq|\.fq.gz|\.fa|\.cram)$", f, flags=re.I
    )
    if bool(match):
        ext = match[1]
        fn = re.sub(rf"\{ext}$", "", f)
    else:
        ext = ""
        fn = f
    # print(fn, ext)
    f = os.path.join("data", "corrected", dataset, fn + ".corrected.fasta")
    return f


def find_input_datasets(wc):
    files = datasets[wc.dataset][wc.type]
    if wc.type == "HQ_herro":
        files = [update_herro_paths(f, wc.dataset) for f in files]
    return files


rule merge_copy_rename_fastq:
    input:
        find_input_datasets,
    output:
        "assembly/input/{dataset}/{dataset}.{type}.fastq.gz",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/merge_copy_rename_fastq.{dataset}.{type}.log",
    shell:
        """
        samtools fastq <(cat {input}) 2>{log}\
        | gzip -c >{output} 2>>{log}
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


rule map_unaligned_bam:
    input:
        bam="data/{path}.bam",
        ref=config["ref"],
    output:
        bam="data/mapped/{path}.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_bam_{path}.log",
    threads: 60
    params:
        preset=lambda wc: "lr:hq" if "duplex" in wc.path else "map-ont",
    shell:
        """
        samtools fastq {input.bam} \
        | minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.ref} - 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

        samtools index {output.bam}
        """


rule bamstats:
    input:
        bam="data/mapped/{path}.bam",
    output:
        hist_l="data/bamstats/{path}/length.hist",
        hist_q="data/bamstats/{path}/quality.hist",
        hist_a="data/bamstats/{path}/accuracy.hist",
        hist_c="data/bamstats/{path}/coverage.hist",
        bamstats="data/bamstats/{path}/bamstats.txt",
        flagstats="data/bamstats/{path}/flagstats.txt",
        basecallers="data/bamstats/{path}/basecallers.txt",
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


rule map_fq:
    input:
        fq="assembly/input/{file}.fastq.gz",
        ref=config["ref"],
    output:
        bam="data/mapped/{file}.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_fq_{file}.log",
    threads: 60
    params:
        preset=lambda wc: "lr:hq" if "duplex" in wc.file else "map-ont",
    shell:
        """
        minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.ref} {input.fq} 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

        samtools index {output.bam}
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

rule sample_to_target_cov:
    input:
        fq="assembly/input/{file}.fastq.gz",
    output:
        fq="assembly/input/{file}.{cov,\d+x}.fastq.gz",
    conda:
        "../env/filtlong.yml"
    log:
        "logs/filtlong_{file}_{cov}.log",
    params:
        min_length=config["min_length"],
        min_mean_q=config["min_mean_q"],
        target_base=lambda wc: str(
            int(wc.cov.replace("x", "")) * config["genome_length"]
        ),
    shell:
        """
        filtlong \
            --min_length {params.min_length} \
            --target_bases {params.target_base} \
            --length_weight 10 \
            --min_mean_q {params.min_mean_q} \
            {input} 2>{log} \
        | gzip -c > {output} 2>>{log}
        """


rule extract_location_data:
    input:
        bam="data/mapped/{file}.bam",
    output:
        fq="assembly/input/{file}.{roi,chr.*}.fastq.gz",
    log:
        "logs/extract_location_data_{file}_{roi}.log",
    shell:
        """
        samtools view -h {input.bam} {wildcards.roi} 2>>{log} \
        | samtools fastq 2>>{log} \
        | gzip -c > {output.fq} 2>>{log}
        """
