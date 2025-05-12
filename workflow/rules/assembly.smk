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

    # Add ULK only if it exists
    if "APK" in datasets[s_d]:
        files["apk"] = f"assembly/input/{s_d}/{s_d}.APK.fastq.gz"

    # Add Trio hapmers if both maternal and paternal are present
    if "HQ_paternal" in datasets[s_d] and "HQ_maternal" in datasets[s_d]:
        if wc.asm in asm_trio.keys():
            files["trio_kmers"] = expand(
                "assembly/input/{dataset}/meryl/{ped}_compress.k30.hapmer.meryl",
                dataset=s_d,
                ped=["maternal", "paternal"],
            )
        elif wc.asm in asm_hifiasm.keys():
            if asm[wc.asm].get("hifiasm", "") == "ont_trio":
                files["trio_kmers"] = expand(
                    "assembly/input/{dataset}/{ped}.yak",
                    dataset=s_d,
                    ped=["maternal", "paternal"],
                )
        else:
            raise ValueError(f"Assembly method {asm} not supported for trio phasing")
    else: files["trio_kmers"] = []
    return files


rule verkko:
    input:
        ul=lambda wc: get_assembly_input(wc).get("ul"),
        hq=lambda wc: get_assembly_input(wc).get("hq"),
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
        "logs/verkko/{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko.txt"
    threads: 92
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


rule copy_verkko_unphased:
    # Input is marked ancient to avoid re-running the rule when verkko scaffolding has finished
    input:
        fa=ancient("assembly/output/verkko/{asm}/assembly.fasta"),
        gfa=ancient("assembly/output/verkko/{asm}/assembly.homopolymer-compressed.gfa"),
        gfa_noseq=ancient(
            "assembly/output/verkko/{asm}/assembly.homopolymer-compressed.noseq.gfa"
        ),
        scfmap=ancient(
            "assembly/output/verkko/{asm}/6-layoutContigs/unitig-popped.layout.scfmap"
        ),
    output:
        fa="assembly/output/verkko_unphased/{asm}/assembly.fasta",
        gfa="assembly/output/verkko_unphased/{asm}/assembly.homopolymer-compressed.gfa",
        gfa_noseq="assembly/output/verkko_unphased/{asm}/assembly.homopolymer-compressed.noseq.gfa",
        scfmap="assembly/output/verkko_unphased/{asm}/assembly.scfmap",
        done="assembly/output/verkko_unphased/{asm}/use_verkko_files.done",
    log:
        "logs/copy_verkko_unphased/{asm}.log",
    shell:
        """
        cp -v {input.fa} {output.fa} > {log} 2>&1
        cp -v {input.gfa} {output.gfa} >> {log} 2>&1
        cp -v {input.gfa_noseq} {output.gfa_noseq} >> {log} 2>&1
        cp -v {input.scfmap} {output.scfmap} >> {log} 2>&1
        touch {output.done}
        """


rule verkko_scaffold:
    input:
        ul=lambda wc: get_assembly_input(wc).get("ul"),
        hq=lambda wc: get_assembly_input(wc).get("hq"),
        porec=lambda wc: get_assembly_input(wc).get("porec"),
        done="assembly/output/verkko_unphased/{asm}/use_verkko_files.done",
    output:
        hp1="assembly/output/verkko/{asm}/assembly.haplotype1.fasta",
        hp2="assembly/output/verkko/{asm}/assembly.haplotype2.fasta",
        colors="assembly/output/verkko/{asm}/assembly.colors.csv",
        scfmap="assembly/output/verkko/{asm}/assembly.scfmap",
        done="assembly/output/verkko/{asm}/create_porec_scaffold.done",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/verkko_scaffold/{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko_scaffold.txt"
    threads: 92
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
        skip_verkko_polish=lambda wc: "--no-correction"
        if asm[wc.asm].get("skip_verkko_polish", True)
        else "",
    shell:
        """
        verkko -d $(dirname {output.hp1}) \
            {params.skip_verkko_polish} \
            --hifi {input.hq} \
            --nano {input.ul} \
            --porec {input.porec} \
            --ahc-run 24 128 6 \
            --fhc-run 24 128 6 \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        touch {output.done}
        """

##############################
## PoreC phasing with GFASE ##
##############################

rule scaffold_create_rename_map:
    input:
        scfmap="assembly/output/verkko_unphased/{asm}/assembly.scfmap",
    output:
        map="assembly/output/gfase/{asm}/contigs.rename.map",
    log:
        "logs/create_rename_map/{asm}.log",
    group:
        "verkko"
    threads: 1
    shell:
        """
        cat {input.scfmap} \
        | grep utig4 \
        | awk '{{print $2"\t"$NF}}' \
        > {output.map} 2>{log}
        """


rule scaffold_rename_fasta:
    input:
        map="assembly/output/gfase/{asm}/contigs.rename.map",
        fa="assembly/output/verkko_unphased/{asm}/assembly.fasta",
    output:
        fa="assembly/output/gfase/{asm}/assembly.fasta",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/scaffold_rename_fasta/{asm}.log",
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
        gfa="assembly/output/verkko_unphased/{asm}/assembly.homopolymer-compressed.gfa",
        fa="assembly/output/gfase/{asm}/assembly.fasta",
    output:
        gfa="assembly/output/gfase/{asm}/assembly.uncompressed.gfa",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/scaffold_uncompress_gfa/{asm}.log",
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
        "logs/scaffold_map_porec/{asm}.log",
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
        "logs/scaffold_gfase/{asm}.log",
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

##################
## Trio phasing ##
##################

# Select correct input reads for kmer counting
def get_trio_input(wc):
    s_d = asm[wc.asm]["dataset"]
    if wc.ped == "child":
        s_hq = "HQ_herro"
    elif wc.ped == "paternal":
        s_hq = "HQ_paternal"
    elif wc.ped == "maternal":
        s_hq = "HQ_maternal"
    return f"assembly/input/{s_d}/{s_d}.{s_hq}.fastq.gz",

# Count homopolymer compressed kmers using meryl on herro corrected reads
rule build_trio_meryldb:
    input:
        get_trio_input
    output:
        meryl=directory("assembly/input/{asm}/meryl/{ped}_compress.k30.meryl"),
    params:
        k=config["K-mer_phasing"],
    conda:
        "../env/merqury.yml"
    log:
        "logs/meryl_builddb/{asm}_{ped}.log",
    threads: 30
    shell:
        """
        meryl count compress\
            k={params.k} \
            threads={threads} \
            {input} \
            output {output.meryl} \
            > {log} 2>&1
        """


rule build_trio_hapmers:
    input:
        expand(
            "assembly/input/{{dataset}}/meryl/{ped}_compress.k30.meryl",
            ped=["maternal", "paternal", "child"],
        ),
    output:
        pat = directory("assembly/input/{dataset}/meryl/paternal_compress.k30.hapmer.meryl"),
        mat = directory("assembly/input/{dataset}/meryl/maternal_compress.k30.hapmer.meryl")
    log:
        "logs/build_trio_hapmers/{dataset}.log",
    threads: 30
    conda:
        "../env/merqury.yml"
    shell:
        """
        # Create an array of input file paths with realpath
        declare -a meryl_dbs=()
        for db in {input}; do
            meryl_dbs+=($(realpath $db))
        done
        
        log=$(realpath {log})

        # Change to output directory
        pushd $(dirname {output[0]})
        
        # Run hapmers script with all input databases
        $MERQURY/trio/hapmers.sh \
            "${{meryl_dbs[@]}}" \
            > $log 2>&1
        
        popd
        """

rule verkko_scaffold_trio:
    input:
        ul=lambda wc: get_assembly_input(wc).get("ul"),
        hq=lambda wc: get_assembly_input(wc).get("hq"),
        kmers=lambda wc: get_assembly_input(wc).get("trio_kmers"),
        done="assembly/output/verkko_unphased/{asm}/use_verkko_files.done",
    output:
        hp1="assembly/output/verkko/{asm}/assembly.haplotype1.fasta",
        hp2="assembly/output/verkko/{asm}/assembly.haplotype2.fasta",
        colors="assembly/output/verkko/{asm}/assembly.colors.csv",
        scfmap="assembly/output/verkko/{asm}/assembly.scfmap",
        done="assembly/output/verkko/{asm}/create_trio_scaffold.done",
    conda:
        "../env/verkko.yml"
    group:
        "verkko"
    log:
        "logs/verkko_scaffold_trio/{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko_scaffold.txt"
    threads: 92
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
        skip_verkko_polish=lambda wc: "--no-correction"
        if asm[wc.asm].get("skip_verkko_polish", False)
        else "",
    shell:
        """
        verkko -d $(dirname {output.hp1}) \
            {params.skip_verkko_polish} \
            --hifi {input.hq} \
            --nano {input.ul} \
            --hap-kmers {input.kmers} trio \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """

ruleorder: verkko_scaffold > verkko_scaffold_trio

#######################
## Hifiasm assembly ##
######################
rule yak_count:
    input:
        get_trio_input
    output:
        yak="assembly/input/{asm}/{ped}.yak"
    conda:
        "../env/yak.yml"
    threads: 24
    benchmark:
        "runtimes/{asm}.yak_count.{ped}.txt"
    log:
        "logs/yak_count/{asm}.{ped}.log"
    shell:
        """
        yak count \
        -t {threads} \
        -k31 \
        -b37 \
        {input} \
        > {output} 2>&1
        """

def get_hifiasm_opts(wildcards, input):
    import traceback
    opts = config["hifiasm_opts"][asm[wildcards.asm].get("hifiasm", "")]
    try:
        opts = opts.format(
            ul = input.ul,
            hq=input.hq,
            maternal = input.kmers[0],
            paternal = input.kmers[1])
    except:
        opts = opts.format(
            ul = input.ul,
            hq=input.hq)
    return opts


rule hifiasm:
    input:
        ul=lambda wc: get_assembly_input(wc).get("ul", []),
        hq=lambda wc: get_assembly_input(wc).get("hq", []),
        kmers=lambda wc: get_assembly_input(wc).get("trio_kmers", []),

    output:
        primary_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.p_ctg.gfa",
        hap1_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.hap1.p_ctg.gfa",
        hap2_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.hap2.p_ctg.gfa",
        primary_gfa_noseq="assembly/output/hifiasm/{asm}/{asm}.bp.p_ctg.noseq.gfa",
        hap1_gfa_noseq="assembly/output/hifiasm/{asm}/{asm}.bp.hap1.p_ctg.noseq.gfa",
        hap2_gfa_noseq="assembly/output/hifiasm/{asm}/{asm}.bp.hap2.p_ctg.noseq.gfa",
    conda:
        "../env/hifiasm.yml"
    log:
        "logs/hifiasm/{asm}.log",
    benchmark:
        "runtimes/{asm}.hifiasm.txt"
    threads: 48
    params:
        hifiasm = config["hifiasm"],
        opts = get_hifiasm_opts
    shell:
        """
        {params.hifiasm} \
            --ul {input.ul} \
            -t {threads} \
            -o $(dirname {output.primary_gfa})/{wildcards.asm} \
            {params.opts} \
        > {log} 2>&1
        """

rule hifiasm_to_fasta:
    input:
        primary_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.p_ctg.gfa",
        hap1_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.hap1.p_ctg.gfa",
        hap2_gfa="assembly/output/hifiasm/{asm}/{asm}.bp.hap2.p_ctg.gfa",
    output:
        fa = "assembly/output/hifiasm/{asm}/assembly.fasta",
        hap1_fa="assembly/output/hifiasm/{asm}/assembly.haplotype1.fasta",
        hap2_fa="assembly/output/hifiasm/{asm}/assembly.haplotype2.fasta",
    log:
        "logs/hifiasm_to_fasta/{asm}.log",
    threads: 1
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.primary_gfa} > {output.fa} 2> {log}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap1_gfa} > {output.hap1_fa} 2>> {log}
        awk '/^S/{{print ">"$2;print $3}}' {input.hap2_gfa} > {output.hap2_fa} 2>> {log}
        """