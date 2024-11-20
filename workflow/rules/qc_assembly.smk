rule bandage:
    input:
        gfa=rules.verkko.output.gfa,
        color=rules.verkko_scaffold.output.colors,
    output:
        svg="assembly/qc/{asm}/bandage_graph.svg",
        png="assembly/qc/{asm}/bandage_graph.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage_{asm}.log",
    threads: 1
    shell:
        """
        Bandage image {input.gfa} {output.svg} --colors {input.color}
        Bandage image {input.gfa} {output.png} --colors {input.color}
        """


rule map_asm_to_ref:
    input:
        fa="assembly/output/{asm}/assembly.{hp}.fasta",
        ref=config["ref"],
    output:
        paf="assembly/qc/{asm}/{hp}.mapped_T2T.paf",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_asm_to_ref.{asm}_{hp}.log",
    threads: 4
    shell:
        """
        minimap2 \
            -cx asm5 \
            -t {threads} \
            {input.ref} {input.fa} \
            > {output.paf} 2> {log}
        """


rule map_cdna_to_asm:
    input:
        asm="assembly/output/{asm}/assembly.{hp}.fasta",
        ref=config["ref_cdna"],
    output:
        paf="assembly/qc/{asm}/cdna_aln.{hp}.paf",
    conda:
        "../env/minimap2.yml"
    threads: 20
    log:
        "logs/qc_asmgene_map_{asm}_{hp}.log",
    shell:
        """
        minimap2 -cxsplice -C5\
            -t {threads} \
            {input.asm} {input.ref} \
            >{output.paf} 2>{log}
        """


rule qc_paftools_stat:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        ref=config["ref"],
    output:
        "assembly/qc/{asm}/qc_paftools_stat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_stat_{asm}_{hp}.log",
    shell:
        """
        paftools.js stat\
            {input.paf} \
            > {output} 2>{log}
        """


rule qc_paftools_asmstat:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        ref=config["ref"],
    output:
        "assembly/qc/{asm}/qc_paftools_asmstat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmstat_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmstat\
            {input.ref}.fai {input.paf} \
            > {output} 2>{log}
        """


rule qc_paftools_asmgene:
    input:
        paf_asm=rules.map_cdna_to_asm.output.paf,
        paf_ref=config["ref_cdna_paf"],
    output:
        "assembly/qc/{asm}/qc_paftools_asmgene.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmgene_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmgene \
            {input.paf_ref} {input.paf_asm} \
            > {output} 2>{log}
        """


rule scaffold_lengths:
    input:
        fa="assembly/output/{asm}/assembly.{hp}.fasta",
        ref=config["ref"],
    output:
        txt="assembly/qc/{asm}/scaffold_lengths.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/scaffold_lengths.{asm}_{hp}.txt",
    shell:
        """
        samtools faidx {input.fa}
        cut -f 1,2 {input.fa}.fai > {output}
        cut -f 1,2 {input.ref}.fai >> {output}
        rm {input.fa}.fai
        """


rule dotplot:
    input:
        len="assembly/qc/{asm}/scaffold_lengths.{hp}.txt",
        paf=rules.map_asm_to_ref.output.paf,
    output:
        "assembly/qc/{asm}/dotplot.{hp}.pdf",
    conda:
        "../env/R.yml"
    threads: 1
    log:
        "logs/dotplot.{asm}_{hp}.log",
    shell:
        """
        Rscript workflow/scripts/minidot.R \
            -i {input.paf} \
            -l {input.len} \
            -o {output} &> {log}
         """


# rule minigraph:
#    input:
#        gfa = rules.verkko.output.gfa,
#        ref = config['ref']
#    output:
#        gfa =
#    conda:
#        "../env/minigraph.yml"
#    shell:
#        """
#        minigraph -cx asm\
#            {input.gfa} {input.ref} \
#            >{output} 2>{log}
#        """

# rule qc_meryl:
#    input:
#
#    output:
#        "assembly_qc/{asm}/kmers.{hp}.meryl"
# rule qc_merqury:
#    input:
#        meryl = expand("assembly_qc/{asm}/kmers.{hp}.meryl", hp = ["haplotype1", "haplotype2"]),
#        fa1 = rules.verkko.output.hp1,
#        fa2 = rules.verkko.output.hp2,
#    output:
#        directory("assembly_qc/{asm}/merqury")
#    conda:
#        "../env/merqury.yml"
#    threads: 1
#    shell:
#    """
#    merqury.sh \
#        {input.meryl} \
#        {input. fa1} {input.fa2} \
#        {output} &> {log}
#    """
