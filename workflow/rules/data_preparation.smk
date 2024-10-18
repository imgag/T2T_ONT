rule map_unaligned_bam:
    input:
        bam = "data/{path}.bam",
        ref = config['ref']
    output:
        bam = "data/mapped/{path}.bam"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_bam_{path}.log"
    threads:
        60
    params:
        preset = lambda wc: "lr:hq"  if "duplex" in wc.path else "map-ont"
    shell:
        """
        samtools fastq {input.bam} \
        | minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.ref} - 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

        samtools index {output}
        """

rule bam_qc:
    input:
        bam="data/mapped/{path}.bam",
        ref=config['ref'],
    output:
        tsv="data/mapped/{path}.qc.tsv.gz",
        json="data/mapped/{path}.qc.json.gz",
    log:
        "logs/bam_qc_{path}.log",
    threads:
        1
    params:
        alfred="bin/alfred"
    shell:
        """
        {params.alfred} qc \
            -r {input.ref} \
            -o {output.tsv} \
            -j {output.json} \
            -i \
            {input.bam} >{log} 2>&1
        """    

rule extract_chr19_testdata:
    input:
        bam_hq = "mapped/published/11_15_22_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.bam",
        bam_ul = "mapped/published/03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup.bam",
        fq_porec = "data/published/SRR27664048.1.1.fastq.gz"
    output:
        fq_hq = "data/published_chr19/published_chr19.HQ.fastq.gz",
        fq_ul = "data/published_chr19/published_chr19.UL.fastq.gz",
        fq_porec = "data/published_chr19/published_chr19.POREC.fastq.gz"
    shell:
        """
        samtools view -h {input.bam_hq} chr19 | samtools fastq -0 {output.fq_hq} -
        samtools view -h {input.bam_ul} chr19 | samtools fastq -0 {output.fq_ul} -
        cp {input.fq_porec} {output.fq_porec}
        """
