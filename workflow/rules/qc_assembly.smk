def get_ref_genome(wc):
    import re

    ref = config["ref"]
    # print(ref)
    match = re.search(r"chr\d+", str(wc.asm))
    # print(match)
    if match:
        # print(match)
        prefix = config["ref"].replace("fasta", "")
        ref = f"{prefix}{match.group(0)}.fasta"
    return ref


def get_assembly_output(wc):
    # Treat undefined isphased as "phased"
    isphased = wc.get("isphased", "phased")

    # Not defined for variant calling

    if isphased == "unphased":
        if wc["tool"] == "verkko":
            return f"assembly/output/verkko_unphased/{wc['asm']}/assembly.fasta"
        else:
            raise ValueError(f"Invalid tool for unphased assembly: {wc['tool']}")

    elif isphased == "phased":
        if wc["tool"] == "verkko":
            if not wc["hp"]:
                return f"assembly/output/verkko/{wc['asm']}/assembly.fasta"
            elif wc["hp"] == "haplotype1":
                return f"assembly/output/verkko/{wc['asm']}/assembly.haplotype1.fasta"
            elif wc["hp"] == "haplotype2":
                return f"assembly/output/verkko/{wc['asm']}/assembly.haplotype2.fasta"
        elif wc["tool"] == "gfase":
            if wc["hp"] == "haplotype1":
                return f"assembly/output/gfase/{wc['asm']}/gfase/phase_0.fasta"
            elif wc["hp"] == "haplotype2":
                return f"assembly/output/gfase/{wc['asm']}/gfase/phase_1.fasta"
        else:
            raise ValueError(f"Invalid tool: {wc['tool']}")



def get_assembly_graph_output(wc):
    if wc["tool"] == "verkko":
        return f"assembly/output/verkko/{wc['asm']}/assembly.homopolymer-compressed.noseq.gfa"
    elif wc["tool"] == "gfase":
        return f"assembly/output/gfase/{wc['asm']}/gfase/chained.gfa"


def get_assembly_graph_colors(wc):
    if wc["tool"] == "verkko":
        return f"assembly/output/verkko/{wc['asm']}/assembly.colors.csv"
    elif wc["tool"] == "gfase":
        return f"assembly/output/gfase/{wc['asm']}/gfase/phases.csv"


rule subsample_ref_genome:
    input:
        fa="data/ref/{ref}.fasta",
    output:
        fa="data/ref/{ref}.{roi,chr.*}.fasta",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/subsample_ref_genome_{ref}_{roi}.log",
    threads: 1
    shell:
        """
        samtools faidx {input.fa} 2>{log}
        samtools faidx {input.fa} {wildcards.roi} > {output} 2>>{log}
        samtools faidx {output.fa} 2>>{log}
        """


rule bandage_unphased:
    input:
        gfa="assembly/output/verkko_unphased/{asm}/assembly.homopolymer-compressed.noseq.gfa",
    output:
        svg="assembly/qc/unphased_verkko/{asm}/bandage_graph.no_colors.svg",
        png="assembly/qc/unphased_verkko/{asm}/bandage_graph.no_colors.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage_unpased_verkko_{asm}.log",
    threads: 1
    shell:
        """
        Bandage image {input.gfa} {output.svg} > {log} 2>&1
        Bandage image {input.gfa} {output.png} >> {log} 2>&1
        """


rule bandage:
    input:
        gfa=get_assembly_graph_output,
        color=rules.verkko_scaffold.output.colors,
    output:
        svg="assembly/qc/phased_{tool}/{asm}/bandage_graph.svg",
        png="assembly/qc/phased_{tool}/{asm}/bandage_graph.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage_{tool}_{asm}.log",
    threads: 1
    shell:
        """
        Bandage image {input.gfa} {output.svg} --colors {input.color} > {log} 2>&1
        Bandage image {input.gfa} {output.png} --colors {input.color} > {log} 2>&1
        """


rule map_asm_to_ref:
    input:
        fa=get_assembly_output,
        ref=get_ref_genome,
    output:
        paf="assembly/qc/{isphased}_{tool}/{asm}/{hp}.mapped_T2T.paf",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_asm_to_ref.{isphased}_{tool}_{asm}_{hp}.log",
    threads: 4
    shell:
        """
        minimap2 \
            -cx asm5 \
            -t {threads} \
            {input.ref} {input.fa} \
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
        "logs/map_cdna_to_ref_{asm}.log",
    shell:
        """
        minimap2 -cxsplice -C5\
            -t {threads} \
            {input.genome} {input.cdna} \
            >{output.paf} 2>{log}
        """


rule map_cdna_to_asm:
    input:
        asm=get_assembly_output,
        ref=config["ref_cdna"],
    output:
        paf="assembly/qc/{isphased}_{tool}/{asm}/cdna_aln.{hp}.paf",
    conda:
        "../env/minimap2.yml"
    threads: 20
    log:
        "logs/qc_asmgene_map_{isphased}_{tool}_{asm}_{hp}.log",
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
        ref=get_ref_genome,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_stat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_stat_{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js stat\
            {input.paf} \
            > {output} 2>{log}
        """


rule qc_paftools_asmstat:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        ref=get_ref_genome,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_asmstat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmstat_{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmstat\
            {input.ref}.fai {input.paf} \
            > {output} 2>{log}
        """


def get_ref_cdna_paf(wc):
    import re

    ref = config["ref"].replace(".fasta", ".cdna.paf")
    # print(ref)
    match = re.search(r"chr\d+", str(wc.asm))
    # print(match)
    if match:
        # print(match)
        prefix = config["ref"].replace("fasta", "")
        ref = f"{prefix}{match.group(0)}.cdna.paf"
    return ref


rule qc_paftools_asmgene:
    input:
        paf_asm=rules.map_cdna_to_asm.output.paf,
        paf_ref=get_ref_cdna_paf,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_asmgene.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmgene_{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmgene \
            {input.paf_ref} {input.paf_asm} \
            > {output} 2>{log}
        """


rule scaffold_lengths:
    input:
        fa=get_assembly_output,
        ref=get_ref_genome,
    output:
        txt="assembly/qc/phased_{tool}/{asm}/scaffold_lengths.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/scaffold_lengths.{tool}_{asm}_{hp}.txt",
    shell:
        """
        samtools faidx {input.fa}
        cut -f 1,2 {input.fa}.fai > {output}
        cut -f 1,2 {input.ref}.fai >> {output}
        rm {input.fa}.fai
        """


rule dotplot:
    input:
        len="assembly/qc/{isphased}_{tool}/{asm}/scaffold_lengths.{hp}.txt",
        paf=rules.map_asm_to_ref.output.paf,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/dotplot.{hp}.pdf",
    conda:
        "../env/R.yml"
    threads: 1
    log:
        "logs/dotplot.{isphased}_{tool}_{asm}_{hp}.log",
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

rule qc_meryl:
    input:
        ref_q100=config["ref_hg002_q100"],
    output:
        meryl=directory("data/ref/hg002_q100_meryl/hg002_q100_k_{k_val}.meryl"),
    params:
        k=config["K-mer"],
    conda:
        "../env/merqury.yml"
    log:
        "logs/meryl_count_ref_q100_{k_val}.log",
    shell:
        """
        meryl count k={params.k} {input.ref_q100} output {output.meryl} > {log} 2>&1
        """


rule qc_merqury_verkko:
    input:
        meryl=f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl',
        pat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype1"}),
        mat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype2"}),
    output:
        out = "assembly/qc/phased_{tool}/{asm}/merqury.qv",
        hap_pat_meryl="assembly/qc/phased_{tool}/{asm}/merqury.haplotype1.qv",
        hap_mat_meryl="assembly/qc/phased_{tool}/{asm}/merqury.haplotype2.qv",
    conda:
        "../env/merqury.yml"
    log:
        "logs/merqury_{tool}_{asm}.log",
    threads: 40
    shell:
        """
        INPUT_MERYL=$(realpath {input.meryl})
        INPUT_PAT_FA=$(realpath {input.pat_fa})
        INPUT_MAT_FA=$(realpath {input.mat_fa})
        OUTPUT_PREFIX=$(dirname $(realpath {output.out}))/merqury
        LOG_FILE=$(realpath {log})
        pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
        cp $INPUT_PAT_FA haplotype1.fa
        cp $INPUT_MAT_FA haplotype2.fa
        export PATH=$PATH:"$CONDA_PREFIX"/share/merqury/eval
        qv.sh \
            $INPUT_MERYL \
            haplotype1.fa \
            haplotype2.fa \
            $OUTPUT_PREFIX \
         >> $LOG_FILE 2>&1
        rm -r *.meryl >> $LOG_FILE 2>&1
        popd >> $LOG_FILE 2>&1
        """


rule qc_merqury_unphased:
    input:
        meryl=f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl',
        fa="assembly/output/gfase/{asm}/assembly.fasta",
    output:
        out="assembly/qc/unphased_verkko/{asm}/merqury.qv",
    conda:
        "../env/merqury.yml"
    log:
        "logs/merqury_unphased_{asm}.log",
    threads: 40
    shell:
        """
        INPUT_MERYL=$(realpath {input.meryl})
        INPUT_PAT_FA=$(realpath {input.fa})
        OUTPUT_PREFIX=$(dirname $(realpath {output.out}))/merqury
        LOG_FILE=$(realpath {log})
        pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
        export PATH=$PATH:"$CONDA_PREFIX"/share/merqury/eval
        qv.sh \
            $INPUT_MERYL \
            $INPUT_PAT_FA \
            $OUTPUT_PREFIX \
         >> $LOG_FILE 2>&1
        rm -r *.meryl >> $LOG_FILE 2>&1
        popd >> $LOG_FILE 2>&1
        """

# temporary meryl files are created during qv.sh process so i put the rm
# export PATH=$PATH:"$CONDA_FREFIX"/share/merqury/eval
# because the qv script

rule find_T2T_contigs:
    input:
        asm = get_assembly_output,
        ref = get_ref_genome,
    output:
        seqinfo = "assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
        alignment = "assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_alignment_T2T.txt",
        motif = "assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif_T2T.txt"
    log:
        "logs/find_T2T_contigs_{isphased}_{tool}_{asm}_{hp}.log",
    threads: 6 
    params:
        T2T_chromosomes = config["T2T_chromosomes"]
    shell:
        """
        export PATH=$PATH:$(dirname {params.T2T_chromosomes})
        {params.T2T_chromosomes} \
            -a {input.asm} \
            -r {input.ref} \
            -o "$(dirname {output.seqinfo})/T2T_contigs.{wildcards.hp}" \
            -m TTAGGG \
            -t {threads} \
            > {log} 2>&1
        """