rule all_plot:
    input:
        expand("results/{sample}/assembly_issues/{sample}.karyogram.gaps.clean.png", sample="T2T00")

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
