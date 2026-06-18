rule all_extended_qc:
    input:        
        # CenMap
        expand("analysis_other/cenmap/{asm}/phased_verkko/cenmap.done", asm = ["T2TP02"]),
        expand("results/{sample}/phased_verkko/ref_alignment/{sample}.synteny.png", sample = ["T2TP02"]),
        expand("results/{sample}/phased_hifiasm/ref_alignment/{sample}.synteny.png", sample = ["T2TP02_hifiasm"]),


# Cenmap
rule cenmap:
    input:
        assembly = lambda wc: get_assembly_output({**wc,'hp': 'both'})["assembly"],
        mod = lambda wc: find_input_datasets(SimpleNamespace(dataset=wc.asm, type="UL"))["files"][0],
        hq = "assembly/input/{asm}/{asm}.HQ_herro.50x.fastq.gz"
    output:
        done = "analysis_other/cenmap/{asm}/{isphased}_{tool}/cenmap.done"
    conda:
        "../env/cenmap.yml"
    log:
        "logs/cenmap/{asm}/cenmap.{isphased}_{tool}.yml"
    threads:
        24
    shell:
        """
        FA=$(realpath {input.assembly})
        HQ=$(realpath {input.hq})
        MOD=$(realpath {input.mod})
        LOG=$(realpath {log})
        WD=$(dirname {output.done})

        pushd $WD >{log}

        cenmap \
            -i $FA \
            -s {wildcards.asm} \
            --hifi $HQ \
            --ont $MOD \
            >$LOG 2>&1
        touch cenmap.done
        """

# Synteny Plot for Genome / Ref alignment
rule install_svbyeye:
  output:
    flag = "bin/.svbyeye_installed"
  log:
    "logs/plot/install_svbyeye.log"
  conda:
    "../env/svbyeye.yml"
  params:
    script = "workflow/scripts/37_install_svbyeye.R"
  shell:
    """
    Rscript {params.script} > {log} 2>&1 && touch {output.flag}
    """

rule synteny_plot_genome:
  input:
    paf_hap1 = "assembly/qc/{isphased}_{tool}/{sample}/haplotype1.mapped_T2T.paf",
    paf_hap2 = "assembly/qc/{isphased}_{tool}/{sample}/haplotype2.mapped_T2T.paf",
    cdna_hap1 = "assembly/qc/{isphased}_{tool}/{sample}/cdna_aln.haplotype1.paf" ,
    cdna_hap2 = "assembly/qc/{isphased}_{tool}/{sample}/cdna_aln.haplotype2.paf" ,
    cdna_ref = config['ref_cdna_paf'],
    flag = "bin/.svbyeye_installed"
  output:
    png = "results/{sample}/{isphased}_{tool}/ref_alignment/{sample}.synteny.png",
    pdf = "results/{sample}/{isphased}_{tool}/ref_alignment/{sample}.synteny.pdf"
  log:
    "logs/plot/{sample}/{isphased}_{tool}/synteny_plot_genome.log"
  threads:
    1
  conda:
    "../env/svbyeye.yml"
  params:
    script = "workflow/scripts/38_synteny_plot.R",
  shell:
    """
    Rscript {params.script} \
        {input.paf_hap1} \
        {input.paf_hap2} \
        {output.png} \
        {output.pdf} \
        {input.cdna_hap1} \
        {input.cdna_hap2} \
        {input.cdna_ref} \
        > {log} 2>&1
    """