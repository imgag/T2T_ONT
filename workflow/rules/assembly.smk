def get_assembly_input(wc):
    s_d = asm[wc.asm]["dataset"]
    s_hq = asm[wc.asm]["HQ_method"]
    s_hq = "HQ_" + s_hq
    s_cov_ul = asm[wc.asm].get("cov_UL", "")
    s_cov_hq = asm[wc.asm].get("cov_HQ", "")
    #print(wc.asm, s_cov_ul, s_cov_hq, s_hq)
    if s_cov_ul != "":
        s_cov_ul = str(s_cov_ul) + "x."
    if s_cov_hq != "":  
        s_cov_hq = str(s_cov_hq) + "x."
    if s_hq == "HQ_combined": 
        s_cov_duplex = asm[wc.asm].get("cov_DUPLEX", "")
        s_cov_herro = asm[wc.asm].get("cov_HERRO", "")
        s_cov_hq = f"{s_cov_duplex}x_{s_cov_herro}x." 
        #print("Combined dataset: ", s_cov_duplex, s_cov_herro)
    s_re = asm[wc.asm].get("region", "")
    if s_re != "":
        s_re = s_re + "."
    files = {
        "ul": f"assembly/input/{s_d}/{s_d}.UL.{s_cov_ul}{s_re}fastq.gz",
        "hq": f"assembly/input/{s_d}/{s_d}.{s_hq}.{s_cov_hq}{s_re}fastq.gz",
        "porec": f"assembly/input/{s_d}/{s_d}.POREC.fastq.gz",
    }
    return files


rule verkko:
    input:
        unpack(get_assembly_input),
    output:
        gfa_noseq="assembly/output/{asm}/assembly.homopolymer-compressed.noseq.gfa",
        gfa="assembly/output/{asm}/assembly.homopolymer-compressed.gfa",
        fa="assembly/output/{asm}/assembly.fasta",
        scfmap = "assembly/output/{asm}/6-layoutContigs/unitig-popped.layout.scfmap"
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/verkko_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko.txt"
    threads: 320
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
    shell:
        """
        verkko -d $(dirname {output.fa}) \
            --hifi {input.hq} \
            --nano {input.ul} \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """

rule verkko_scaffold:
    input:
        unpack(get_assembly_input),
    output:
        hp1="assembly/output/{asm}/assembly.haplotype1.fasta",
        hp2="assembly/output/{asm}/assembly.haplotype2.fasta",
        colors="assembly/output/{asm}/assembly.colors.csv",
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/verkko_scaffold_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko_scaffold.txt"
    threads: 120
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
    shell:
        """
        verkko -d $(dirname {output.hp1}) \
            --hifi {input.hq} \
            --nano {input.ul} \
            --porec {input.porec} \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """

rule scaffold_create_rename_map:
    input:
        scfmap = rules.verkko.output.scfmap
    output:
        "assembly/scaffold/{asm}/contigs.rename.map"
    log:
        "logs/create_rename_map_{asm}.log"
    group:
        "verrko"
    threads:
        1
    shell:
        """
        cat {input.scfmap} \
        | grep utig4 \
        | awk '{{print $2"\t"$NF}}' \
        > {output} 2>{log}
        """

rule scaffold_rename_fasta:
    input:
        map = "assembly/scaffold/{asm}/contigs.rename.map",
        fa = "assembly/output/{asm}/assembly.fasta"
    output:
        fa = "assembly/scaffold/{asm}/assembly.fasta"
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/scaffold_rename_fasta_{asm}.log"
    threads:
        1
    shell:
        """
        $CONDA_PREFIX/lib/verkko/scripts/fasta_combine.py rename\
            {output.fa} \
            {input.map} \
            {input.fa} \
            >{log} 2>&1
        """

rule scaffold_uncompress_gfa:
    input:
        gfa = "assembly/output/{asm}/assembly.homopolymer-compressed.gfa",
        fa = "assembly/scaffold/{asm}/assembly.fasta"
    output:
        gfa = "assembly/scaffold/{asm}/assembly.uncompressed.gfa"
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/scaffold_uncompress_gfa_{asm}.log"
    threads:
        4
    shell:
        """
        $CONDA_PREFIX/lib/verkko/bin/alignGFA \
            -V -e 0.30 \
            -gfa \
            -i {input.gfa} \
            -T {input.fa} 0 \
            -t {threads} \
            -o {output.gfa} \
            >{log} 2>&1
        """

rule scaffold_map_porec:
    input:
        unpack(get_assembly_input),
        asm = "assembly/scaffold/{asm}/assembly.fasta"
    output:
        bam = "assembly/scaffold/{asm}/asm_porec.bam"
    conda:
        "../env/minimap2.yml"
    log:
        "logs/scaffold_map_porec_{asm}.log"
    threads:
        60
    shell:
        """
        minimap2 \
            -a \
            -x map-ont \
            -k 17 \
            -t 56 \
            -K 10g \
            -I 8g \
            {input.asm} \
            {input.porec} \
            2>{log} \
        | samtools view -bh -@ 8 -q 1 - \
            > {output} 2>>{log}
        """

rule scaffold_gfase:
    input:
        bam_porec = "assembly/scaffold/{asm}/asm_porec.bam",
        gfa = "assembly/scaffold/{asm}/assembly.uncompressed.gfa"
    output:
        directory("assembly/scaffold/{asm}/gfase")
    params:
        gfase = "bin/GFAse/build/phase_contacts_with_monte_carlo"
    log:
        "logs/scaffold_gfase_{asm}.log"
    threads:
        62
    shell:
        """
        {params.gfase} \
            -i {input.bam_porec} \
            -g {input.gfa} \
            -o {output} \
            --use_homology \
            --skip_unzip \
            -m 3 \
            -t {threads} >{log} 2>&1
        """