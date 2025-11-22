from types import SimpleNamespace

# rule map_unaligned_bam:
#     input:
#         bam="data/{path}.bam",
#         ref=config["ref"],
#     output:
#         bam="data/mapped/{path}.ref.bam",
#     conda:
#         "../env/minimap2.yml"
#     log:
#         "logs/map_bam/{path}.log",
#     threads: 60
#     params:
#         preset=lambda wc: "lr:hq" if "duplex" in wc.path else "map-ont",
#     shell:
#         """
#         samtools fastq {input.bam} \
#         | minimap2 --MD -ax {params.preset} --eqx \
#             -t {threads} \
#             {input.ref} - 2>{log} \
#         | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

#         samtools index {output.bam}
#         """

# Map fq from input folder or correction to ref genome for QC, coverage and bamstats
rule map_fq_ref:
    input:
        fq="assembly/input/{file}.fastq.gz",
        ref= config["ref"],
    output:
        bam="data/mapped/{file}.ref.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_fq/{file}.log",
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

rule create_fq_for_map_with_mod:
    input:
        bam = lambda wc: find_input_datasets(SimpleNamespace(dataset=wc.asm, type="UL"))["files"],
        ref = config["ref"],
    output:
        fq = temp("data/mapped/{asm}/{asm}.UL.100G.mod.fastq")
    log:
        "logs/mapping/{asm}/create_fq_for_map_with_mod.log"
    conda:
        "../env/samtools.yml"
    threads:    
        1
    shell:
        """
        samtools fastq -TMM,ML --reference {input.ref} <(samtools cat {input.bam}) \
        | head -c 1000GB \
        > {output.fq}.tmp 2>>{log}
        
        lines=$(wc -l < {output.fq}.tmp)
        complete_records=$((lines / 4 * 4))
        head -n $complete_records {output.fq}.tmp > {output.fq}
        rm {output.fq}.tmp
        """

# Map UL (including modifications to diploid assembly)
rule map_ul_to_asm:
    input:
        fa = lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
            if wc.ref == "asm" else config["ref"],
        fq = "data/mapped/{asm}/{asm}.UL.100G.mod.fastq"
    output:
        bam = "data/mapped/{asm}/{asm}.UL.100G.{ref}.bam"
    wildcard_constraints:
        ref="[^.]+"
    log:
        "logs/mapping/{asm}/map_hq_to_asm.{ref}.log"
    conda:
        "../env/minimap2.yml"
    threads:
        40
    shell:
        # Convert to FASTQ and keep modifications
        # Remove unmapped, non primary, secondary alignments
        """
        minimap2 --MD -ax lr:hq --eqx -y \
            -t {threads} \
            {input.fa} {input.fq} 2>>{log} \
        | samtools view -h -F 2308 \
        | samtools sort -m 8G -@ 8 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """

# Map HQ to diploid assembly, for Nucflag and Flagger
rule map_hq_to_asm:
    input:
        assembly = lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"],
        hq = lambda wc: get_assembly_input(wc).get("hq")
    output:
        bam = "data/mapped/{asm}/{asm}.HQ.asm.bam"
    log:
        "logs/mapping/{asm}/map_hq_to_asm.log"
    conda:
        "../env/minimap2.yml"
    threads:
        40
    params:
        preset="lr:hqae",
    shell:
        """
        minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.assembly} {input.hq} 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1
        samtools index {output.bam}
        """


rule map_asm_to_ref:
    input:
        unpack(lambda wc: {k: ancient(v) for k, v in get_assembly_output(wc).items()}),
        ref=get_ref_genome,
        phased_out=lambda wc: (
            f"assembly/output/verkko/{wc.asm}/assembly.colors.csv"
            if wc.isphased == "phased" and wc.tool == "verkko"
            else []
        ),
    output:
        paf="assembly/qc/{isphased}_{tool}/{asm}/{hp}.mapped_T2T.paf",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_asm_to_ref/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 4
    shell:
        """
        minimap2 \
            -cx asm5 \
            -t {threads} \
            {input.ref} {input.assembly} \
            > {output.paf} 2> {log}
        """


rule map_cdna_to_ref:
    input:
        genome=get_ref_genome,
        cdna=config["ref_cdna"],
    output:
        paf="data/ref/{asm}.cdna.paf",
    conda:
        "../env/minimap2.yml"
    threads: 20
    log:
        "logs/map_cdna_to_ref/{asm}.log",
    shell:
        """
        minimap2 -cxsplice -C5\
            -t {threads} \
            {input.genome} {input.cdna} \
            >{output.paf} 2>{log}
        """


rule map_cdna_to_asm:
    input:
        unpack(get_assembly_output),
        ref=config["ref_cdna"],
    output:
        paf="assembly/qc/{isphased}_{tool}/{asm}/cdna_aln.{hp}.paf",
    conda:
        "../env/minimap2.yml"
    threads: 20
    log:
        "logs/qc_asmgene_map/{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        minimap2 -cxsplice -C5\
            -t {threads} \
            {input.assembly} {input.ref} \
            >{output.paf} 2>{log}
        """
