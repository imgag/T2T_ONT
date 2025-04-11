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
            if not wc["hp"] or wc["hp"] == "both":
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


def get_assembly_colors(wc):
    if wc["isphased"] == "phased":
        if wc["tool"] == "verkko":
            return f"assembly/output/verkko/{wc['asm']}/assembly.colors.csv"
        if wc["tool"] == "gfase":
            return f"assembly/output/gfase/{wc['asm']}/gfase/phases.csv"
    elif wc["isphased"] == "unphased":
        return f"assembly/qc/unphased_{wc['tool']}/{wc['asm']}/colors.tsv"
    else:
        raise ValueError(f"Invalid  phasing value: {wc['isphased']}")


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


rule create_colors:
    input:
        paf="assembly/qc/{isphased}_{tool}/{asm}/both.mapped_T2T.paf",
        colors_phasing=lambda wc: f"assembly/output/verkko/{wc.asm}/assembly.colors.csv"
        if wc.isphased == "phased" and wc.tool == "verkko"
        else [],
    output:
        csv="assembly/qc/{isphased}_{tool}/{asm}/colors.tsv",
    log:
        "logs/create_colors_{isphased}_{tool}_{asm}.log",
    params:
        colours_phasing=lambda wc: f"-c assembly/output/verkko/{wc.asm}/assembly.colors.csv"
        if wc.isphased == "phased" and wc.tool == "verkko"
        else "",
    shell:
        """
        python workflow/scripts/11_extract_colors.py \
            -i {input.paf} {params.colours_phasing} \
            -o {output.csv} \
            >{log} 2>&1
        """


rule process_graph:
    input:
        gfa=get_assembly_graph_output,
        scfmap="assembly/output/verkko/{asm}/assembly.scfmap",
        color="assembly/qc/{isphased}_{tool}/{asm}/colors.tsv",
    output:
        gfa="assembly/qc/{isphased}_{tool}/{asm}/assembly_graph.gfa",
    log:
        "logs/process_graph_{isphased}_{tool}_{asm}.log",
    shell:
        """
        python workflow/scripts/12_process_gfa.py \
            --gfa {input.gfa} \
            --scfmap {input.scfmap} \
            --colors {input.color} \
            --output {output.gfa} \
            > {log} 2>&1
        """


rule process_graph_phased:
    input:
        gfa=get_assembly_graph_output,
        scfmap="assembly/output/verkko_unphased/{asm}/assembly.scfmap",
        color="assembly/qc/{isphased}_{tool}/{asm}/colors.tsv",


rule bandage:
    input:
        gfa="assembly/qc/{isphased}_{tool}/{asm}/assembly_graph.gfa",
        color=get_assembly_colors,
    output:
        svg="assembly/qc/{isphased}_{tool}/{asm}/bandage_graph.svg",
        png="assembly/qc/{isphased}_{tool}/{asm}/bandage_graph.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage_{isphased}_{tool}_{asm}.log",
    threads: 1
    shell:
        """
        Bandage image {input.gfa} {output.svg} --colors {input.color} > {log} 2>&1
        Bandage image {input.gfa} {output.png} --colors {input.color} > {log} 2>&1
        """


rule map_asm_to_ref:
    input:
        fa=ancient(get_assembly_output),
        ref=get_ref_genome,
        phased_out=lambda wc: f"assembly/output/verkko/{wc.asm}/assembly.colors.csv"
        if wc.isphased == "phased" and wc.tool == "verkko"
        else [],
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


# Select correct Meryl ref for sample. Two options:
# 1) HG002 for published set, fallback 2) Illumina shortread data for TUE_02 sample
def get_meryl_ref(wc):
    if "published" in str(wc.asm):
        return f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl'
    elif "TUE_02" in str(wc.asm):
        return "analysis_other/merqury_shortread/DX203429_02.meryl"
    else:
        print(f"WARNING: Unknown Meryl Ref for {wc.asm}, using HG002")
        return f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl'


rule qc_merqury_phased:
    input:
        meryl=get_meryl_ref,
        pat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype1"}),
        mat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype2"}),
    output:
        out="assembly/qc/phased_{tool}/{asm}/merqury.qv",
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


def get_merqury_input(wc):
    if wc["isphased"] == "unphased":
        return f"assembly/output/verkko_unphased/{wc['asm']}/assembly{wc['polished']}.fasta"
    elif wc["isphased"] == "phased":
        return f"assembly/output/verkko/{wc['asm']}/assembly{wc['polished']}.fasta"


rule qc_merqury_unphased:
    input:
        meryl=get_meryl_ref,
        fa=get_merqury_input,
    output:
        out="assembly/qc/{isphased}_verkko/{asm}/merqury{polished}.qv",
    wildcard_constraints:
        polished=".*",
    conda:
        "../env/merqury.yml"
    log:
        "logs/merqury_{isphased}_verkko_{asm}{polished}.log",
    threads: 40
    shell:
        """
        INPUT_MERYL=$(realpath {input.meryl})
        INPUT_PAT_FA=$(realpath {input.fa})
        OUTPUT_PREFIX=$(dirname $(realpath {output.out}))/merqury{wildcards.polished}
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


ruleorder: qc_merqury_phased > qc_merqury_unphased


# temporary meryl files are created during qv.sh process so i put the rm
# export PATH=$PATH:"$CONDA_FREFIX"/share/merqury/eval
# because the qv script


rule find_T2T_contigs:
    input:
        asm=get_assembly_output,
        ref=get_ref_genome,
    output:
        seqinfo="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
        alignment="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_alignment_T2T.txt",
        motif="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif_T2T.txt",
    log:
        "logs/find_T2T_contigs_{isphased}_{tool}_{asm}_{hp}.log",
    threads: 6
    params:
        T2T_chromosomes=config["T2T_chromosomes"],
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
