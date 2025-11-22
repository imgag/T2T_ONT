rule all_methylation:
    input:
        expand("analysis_other/mod/{sample}/{sample}.{ref}.modkit_summary.tsv", sample = finished_samples, ref = ["asm", "ref"]),
        expand("analysis_other/mod/{sample}/{sample}.{ref}.methyl.bed", sample = finished_samples, ref = ["asm", "ref"]),
        expand("analysis_other/mod/{sample}/{sample}.{ref}.bw", sample = finished_samples, ref = ["asm", "ref"]),
        expand("analysis_other/mod/{sample}/{sample}.ref.hp{hp}.methyl.bed", sample = finished_samples, hp = ["1", "2", "0"]),
        expand("analysis_other/mod/{sample}/{sample}.ref.hp{hp}.bw", sample = finished_samples, hp = ["1", "2", "0"])

rule modkit_summary:
    input:
        "data/mapped/{sample}/{sample}.UL.100G.{ref}.bam"
    output:
        "analysis_other/mod/{sample}/{sample}.{ref}.modkit_summary.tsv"
    log:
        "logs/methylation/{sample}.modkit_summary.{ref}.log"
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


rule longphase_haplotag:
    input:
        bam = "data/mapped/{sample}/{sample}.UL.100G.ref.bam",
        ref = config['ref'],
        vcf = "assembly/variants/{sample}/phased_verkko/small_variants.dip.vcf.gz"
    output:
        "data/mapped/{sample}/{sample}.UL.100G.haplotagged.ref.bam"
    params:
        longphase = config['longphase']
    threads:
        4
    shell:
        """
        {params.longphase} haplotag \
            --snp-file {input.vcf} \
            --bam-file {input.bam} \
            --reference {input.ref} \
            --threads {threads} \
            --out-prefix data/mapped/{wildcards.sample}/{wildcards.sample}.UL.70x.ref.haplotagged \
            >{log} 2>&1

        samtools index {output} \
            >{log} 2>&1
        """

rule modkit_pileup_complete:
    input:
        bam = "data/mapped/{sample}/{sample}.UL.100G.{ref}.bam",
        ref = config['ref']
    output:
        "analysis_other/mod/{sample}/{sample}.{ref}.methyl.bed"
    log:
        "logs/methylation/{sample}.{ref}.modkit_pileup_complete.log"
    params:
        modkit = config['modkit']
    resources:
        mem_gb = 150
    threads:
        8
    shell:
        """
        {params.modkit} pileup \
            --threads {threads} \
            --ref {input.ref} \
            {input.bam} {output} \
            >{log} 2>&1
        """

rule extract_haplotype_reads:
    input:
        bam = "data/mapped/{sample}/{sample}.UL.100G.haplotagged.ref.bam"
    output:
        bam = temp("data/mapped/{sample}/{sample}.UL.100G.ref.hp{hp}.bam"),
        bai = temp("data/mapped/{sample}/{sample}.UL.100G.ref.hp{hp}.bam.bai")
    log:
        "logs/methylation/{sample}.extract_hp{hp}.log"
    threads:
        4
    wildcard_constraints:
        haplotype = "HP1|HP2|HP0"
    shell:
        """
        samtools view \
            -@ {threads} \
            -h -d HP:{wildcards.hp} \
            {input.bam} \
            -b -o {output.bam} \
            >{log} 2>&1
        
        samtools index {output.bam} \
            >>{log} 2>&1
        """

rule modkit_pileup_haplotype:
    input:
        bam = "data/mapped/{sample}/{sample}.UL.100G.ref.hp{hp}.bam",
        ref = config['ref']
    output:
        "analysis_other/mod/{sample}/{sample}.hp{hp}.methyl.bed"
    log:
        "logs/methylation/{sample}.modkit_pileup_hp{hp}.log"
    params:
        modkit = config['modkit']
    threads:
        8
    wildcard_constraints:
        haplotype = "0|1|2"
    shell:
        """
        {params.modkit} pileup \
            --threads {threads} \
            --ref {input.ref} \
            --log-filepath {log} \
            {input.bam} {output} \
            2>>{log}
        """

rule bedmethyl_to_bigwig:
    input:
        bedmethyl = "analysis_other/mod/{sample}/{file}.methyl.bed",
        chrom_sizes = config['ref'] + ".fai"
    output:
        "analysis_other/mod/{sample}/{file}.bw"
    log:
        "logs/methylation/{sample}.bedmethyl_to_bigwig.{file}.log"
    conda:
        "../env/ucsc-binaries.yml"
    threads:
        1
    shell:
        """
        bedGraphToBigWig \
            <(awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $4/($4+$5)}}' {input.bedmethyl} ) \
            {input.chrom_sizes} \
            {output} \
        >{log} 2>&1
        """

