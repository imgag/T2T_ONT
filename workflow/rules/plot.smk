rule all_plot:
    input:
        expand("results/{sample}/assembly_issues/{sample}.karyogram.gaps.clean.png", sample="T2T00"),
        expand("results/{sample}/ref_alignment/{sample}.synteny.png", sample=finished_samples)

# Karyogram plotting rule for Snakemake
rule plot_karyogram:
    input:
        fai=config['ref']+".fai",
        seqinfo="assembly/qc/phased_verkko/{sample}/T2T_contigs.both.seqinfo.txt",
        colors="assembly/qc/phased_verkko/{sample}/colors.tsv",
        gap_stats="assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed",
        nucflag="analysis_other/nucflag/{sample}/nucflag_misasm.bed",
        rm_satellites="analysis_other/repeatmasker/{sample}/rm_summary/{sample}_filtered_satellites.bed",
        telo="assembly/qc/phased_verkko/{sample}/T2T_contigs.both_telomeric_repeat_windows.csv"
    output:
        gaps_all="results/{sample}/assembly_issues/{sample}.karyogram.gaps.with_unassigned.png",
        gaps_clean="results/{sample}/assembly_issues/{sample}.karyogram.gaps.clean.png",
        nucflag_clean="results/{sample}/assembly_issues/{sample}.karyogram.nucflag.clean.png",
        gap_lengths="results/{sample}/assembly_issues/{sample}.gap_lengths.png",
        nucflag_types="results/{sample}/assembly_issues/{sample}.nucflag_types.png"
    params:
        script="workflow/scripts/36_plot_karyograms.R",
        output_prefix="results/{sample}/assembly_issues/{sample}"
    conda:
        "../env/R.yml"
    threads: 1
    log:
        "logs/plot/{sample}/plot_karyograms.log"
    shell:
        """
        Rscript {params.script} \
            {input.fai} \
            {input.seqinfo} \
            {input.colors} \
            {input.gap_stats} \
            {input.nucflag} \
            {input.rm_satellites} \
            {input.telo} \
            {params.output_prefix} \
            > {log} 2>&1
        """


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
    paf_hap1 = "assembly/qc/phased_verkko/{sample}/haplotype1.mapped_T2T.paf",
    paf_hap2 = "assembly/qc/phased_verkko/{sample}/haplotype2.mapped_T2T.paf",
    cdna_hap1 = "assembly/qc/phased_verkko/{sample}/cdna_aln.haplotype1.paf" ,
    cdna_hap2 = "assembly/qc/phased_verkko/{sample}/cdna_aln.haplotype2.paf" ,
    cdna_ref = config['ref_cdna_paf'],
    flag = "bin/.svbyeye_installed"
  output:
    png = "results/{sample}/ref_alignment/{sample}.synteny.png",
    pdf = "results/{sample}/ref_alignment/{sample}.synteny.pdf"
  log:
    "logs/plot/{sample}/synteny_plot_genome.log"
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