chromosome_colors = {
    # Main chromosomes - using colorblind-friendly base colors with light/dark variants
    'chr1':  {'haplotype1': '#E69F00', 'haplotype2': '#FFB319'},
    'chr2':  {'haplotype1': '#56B4E9', 'haplotype2': '#7BC6EE'},
    'chr3':  {'haplotype1': '#009E73', 'haplotype2': '#00BF8C'},
    'chr4':  {'haplotype1': '#CC79A7', 'haplotype2': '#D694B8'},
    'chr5':  {'haplotype1': '#0072B2', 'haplotype2': '#0089D9'},
    'chr6':  {'haplotype1': '#D55E00', 'haplotype2': '#FF7400'},
    'chr7':  {'haplotype1': '#666666', 'haplotype2': '#999999'},
    
    # Recycling colors with different shades
    'chr8':  {'haplotype1': '#E6AB23', 'haplotype2': '#FFCD66'},
    'chr9':  {'haplotype1': '#5699E9', 'haplotype2': '#89B9F0'},
    'chr10': {'haplotype1': '#00A880', 'haplotype2': '#33BF99'},
    'chr11': {'haplotype1': '#CC8DB3', 'haplotype2': '#D6A8C4'},
    'chr12': {'haplotype1': '#1A7FBA', 'haplotype2': '#4D9ECC'},
    'chr13': {'haplotype1': '#D57533', 'haplotype2': '#FF9B66'},
    'chr14': {'haplotype1': '#737373', 'haplotype2': '#A6A6A6'},
    'chr15': {'haplotype1': '#E6B847', 'haplotype2': '#FFD480'},
    'chr16': {'haplotype1': '#567DE9', 'haplotype2': '#89A3F0'},
    'chr17': {'haplotype1': '#00B28C', 'haplotype2': '#33C6A6'},
    'chr18': {'haplotype1': '#CCA0BF', 'haplotype2': '#D6BBD0'},
    'chr19': {'haplotype1': '#338CC2', 'haplotype2': '#66ACD9'},
    'chr20': {'haplotype1': '#D58C66', 'haplotype2': '#FFB499'},
    'chr21': {'haplotype1': '#808080', 'haplotype2': '#B3B3B3'},
    
    # Sex chromosomes
    'chrX':  {'haplotype1': '#9467BD', 'haplotype2': '#B189D6'},
    'chrY':  {'haplotype1': '#8C564B', 'haplotype2': '#A67C73'}
}

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


def get_assembly_colors(wc):
    if wc["isphased"] == "phased":
        if wc["tool"] == "verkko":
            return f"assembly/output/verkko/{wc['asm']}/assembly.colors.csv"
        if wc["tool"] == "gfase":
            return f"assembly/output/gfase/{wc['asm']}/gfase/phases.csv"
    elif wc["isphased"] == "unphased":
        return f"assembly/qc/unphased_{wc['tool']}/{wc['asm']}/colors.csv"
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
        paf="assembly/qc/{isphased}_{tool}/{asm}/{isphased}.mapped_T2T.paf"
    output: 
        csv = "assembly/qc/{isphased}_{tool}/{asm}/colors.csv"
    run:
        contig_colors = []
        hp = "haplotype1"
        with open(input.paf, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                contig_id = fields[0]
                chrom = fields[5]
                color = chromosome_colors[contig_id][hp]
                new_line = f"{contig_id}\t{color}"
                contig_colors.append(new_line)
        with open(output.csv, "w") as f:
            f.write("contig\tcolor\n")
            f.write("\n".join(contig_colors))

rule bandage:
    input:
        gfa=get_assembly_graph_output,
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