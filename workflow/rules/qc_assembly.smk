rule bandage:
    input:
        gfa = rules.verkko.output.gfa
    output:
        svg = "assembly_qc/{dataset}/bandage_graph.svg",
        png = "assembly_qc/{dataset}/bandage_graph.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage_{dataset}.log"
    threads:
        1
    shell:
        """
        Bandage image {input.gfa} {output.svg}
        Bandage image {input.gfa} {output.png}
        """

rule map_asm_to_ref:
    input:
        fa = "assembly/{dataset}/assembly.{hp}.fasta",
        ref = config['ref']
    output:
        paf = "assembly_qc/{dataset}/{hp}.mapped_T2T.paf"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_asm_to_ref.{dataset}_{hp}.log"
    threads:
        4
    shell:
        """
        minimap2 \
            -cx asm5 \
            -t {threads} \
            {input.ref} {input.fa} \
            > {output.paf} 2> {log}
        """

rule qc_paftools:
    input:
        paf = rules.map_asm_to_ref.output.paf,
        ref = config['ref']
    output:
        "assembly_qc/{dataset}/qc_paftools.{hp}.txt"
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_{dataset}_{hp}.log"
    shell:
        """
        paftools.js asmstat {input.ref}.fai {input.paf} > {output}
        """

rule scaffold_lengths:
    input:
        fa = "assembly/{dataset}/assembly.{hp}.fasta",
        ref = config['ref']
    output:
        txt = "assembly_qc/{dataset}/scaffold_lengths.{hp}.txt"
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/scaffold_lengths.{dataset}_{hp}.txt"
    shell:
        """
        samtools faidx {input.fa}
        cut -f 1,2 {input.fa}.fai > {output}
        cut -f 1,2 {input.ref}.fai >> {output}
        rm {input}.fai
        """
    
rule dotplot:
    input:
        len = "assembly_qc/{dataset}/scaffold_lengths.{hp}.txt",
        paf = rules.map_asm_to_ref.output.paf
    output:
        "assembly_qc/{dataset}/dotplot.{hp}.pdf"
    conda:
        "../env/R.yml"
    threads:
        1
    log:
        "logs/dotplot.{dataset}_{hp}.log"
    shell:
        """
        Rscript workflow/scripts/minidot.R \
            -i {input.paf} \
            -l {input.len} \
            -o {output} &> {log}
         """

#rule qc_meryl:
#    input:
#        
#    output:
#        "assembly_qc/{dataset}/kmers.{hp}.meryl"

#rule qc_merqury:
#    input:
#        meryl = expand("assembly_qc/{dataset}/kmers.{hp}.meryl", hp = ["haplotype1", "haplotype2"]),
#        fa1 = rules.verkko.output.hp1,
#        fa2 = rules.verkko.output.hp2,
#    output:
#        directory("assembly_qc/{dataset}/merqury")
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