def get_assembly_input(wc):
    # Get the dataset, HQ method, and coverage from the assembly config
    s_d = asm[wc.asm]["dataset"]
    s_hq = asm[wc.asm]["HQ_method"]
    s_hq = "HQ_" + s_hq
    s_cov_ul = asm[wc.asm].get("cov_UL", "")
    s_cov_hq = asm[wc.asm].get("cov_HQ", "")
    # print(wc.asm, s_cov_ul, s_cov_hq, s_hq)
    # Handle coverage and use full dataset if no coverage is specified
    if s_cov_ul != "":
        s_cov_ul = str(s_cov_ul) + "x."
    if s_cov_hq != "":
        s_cov_hq = str(s_cov_hq) + "x."
    if s_hq == "HQ_combined":
        s_cov_duplex = asm[wc.asm].get("cov_DUPLEX", "")
        s_cov_herro = asm[wc.asm].get("cov_HERRO", "")
        s_cov_hq = f"{s_cov_duplex}x_{s_cov_herro}x."
        # print("Combined dataset: ", s_cov_duplex, s_cov_herro)
    s_re = asm[wc.asm].get("region", "")
    if s_re != "":
        s_re = s_re + "."
    files = {
        "ul": f"assembly/input/{s_d}/{s_d}.UL.{s_cov_ul}{s_re}fastq.gz",
        "hq": f"assembly/input/{s_d}/{s_d}.{s_hq}.{s_cov_hq}{s_re}fastq.gz",
    }
    # Add POREC only if it exists
    if "POREC" in datasets[s_d]:
        files["porec"] = f"assembly/input/{s_d}/{s_d}.POREC.fastq.gz"
    return files


rule verkko:
    input:
        unpack(get_assembly_input),
    output:
        gfa_noseq="assembly/output/verkko/{asm}/assembly.homopolymer-compressed.noseq.gfa",
        gfa="assembly/output/verkko/{asm}/assembly.homopolymer-compressed.gfa",
        fa="assembly/output/verkko/{asm}/assembly.fasta",
        scfmap="assembly/output/verkko/{asm}/6-layoutContigs/unitig-popped.layout.scfmap",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/verkko_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko.txt"
    threads: 30
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


#rule copy_verrko_unphased:
#    input: 

rule verkko_scaffold:
    input:
        unpack(get_assembly_input),
        done="assembly/output/gfase/{asm}/use_verkko_files.done",
    output:
        hp1="assembly/output/verkko/{asm}/assembly.haplotype1.fasta",
        hp2="assembly/output/verkko/{asm}/assembly.haplotype2.fasta",
        colors="assembly/output/verkko/{asm}/assembly.colors.csv",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/verkko_scaffold_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko_scaffold.txt"
    threads: 30
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
        scfmap=rules.verkko.output.scfmap,
    output:
        "assembly/output/gfase/{asm}/contigs.rename.map",
    log:
        "logs/create_rename_map_{asm}.log",
    group:
        "verkko"
    threads: 1
    shell:
        """
        cat {input.scfmap} \
        | grep utig4 \
        | awk '{{print $2"\t"$NF}}' \
        > {output} 2>{log}
        """


rule scaffold_rename_fasta:
    input:
        map="assembly/output/gfase/{asm}/contigs.rename.map",
        fa=ancient("assembly/output/verkko/{asm}/assembly.fasta"),
    output:
        fa="assembly/output/gfase/{asm}/assembly.fasta",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/scaffold_rename_fasta_{asm}.log",
    threads: 1
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
        gfa=ancient("assembly/output/verkko/{asm}/assembly.homopolymer-compressed.gfa"),
        fa="assembly/output/gfase/{asm}/assembly.fasta",
    output:
        gfa="assembly/output/gfase/{asm}/assembly.uncompressed.gfa",
        done="assembly/output/gfase/{asm}/use_verkko_files.done",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/scaffold_uncompress_gfa_{asm}.log",
    threads: 4
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
        
        touch {output.done}
        """


rule scaffold_map_porec:
    input:
        unpack(get_assembly_input),
        asm="assembly/output/gfase/{asm}/assembly.fasta",
    output:
        bam="assembly/output/gfase/{asm}/asm_porec.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/scaffold_map_porec_{asm}.log",
    threads: 40
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
        bam_porec="assembly/output/gfase/{asm}/asm_porec.bam",
        gfa="assembly/output/gfase/{asm}/assembly.uncompressed.gfa",
    output:
        asm_hp1="assembly/output/gfase/{asm}/gfase/phase_0.fasta",
        asm_hp2="assembly/output/gfase/{asm}/gfase/phase_1.fasta",
        gfa="assembly/output/gfase/{asm}/gfase/chained.gfa",
    params:
        gfase="bin/GFAse/build/phase_contacts_with_monte_carlo",
    log:
        "logs/scaffold_gfase_{asm}.log",
    threads: 12
    shell:
        """
        outdir=$(dirname {output.gfa})
        rm -rf $outdir
        {params.gfase} \
            -i {input.bam_porec} \
            -g {input.gfa} \
            -o  $outdir \
            --use_homology \
            --skip_unzip \
            -m 3 \
            -t {threads} >{log} 2>&1
        """
