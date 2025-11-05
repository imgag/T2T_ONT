# Karyogram plotting rule for Snakemake
# Generates chromosome QC plots showing assembly quality, gaps, repeats, and misassemblies

rule plot_karyogram:
    input:
        fai="data/ref/{reference}.fasta.fai",
        seqinfo="assembly/qc/{assembler}/{sample}/T2T_contigs.unphased.seqinfo.txt",
        colors="assembly/qc/{assembler}/{sample}/colors.tsv",
        gap_stats="assembly/qc/{assembler}/{sample}/gap_stats.both.n_regions.bed",
        nucflag="analysis_other/nucflag/{sample}.nucflag_misasm.bed",
        rm_satellites="analysis_other/repeatmasker/{sample}/rm_summary_filtered_satellites.bed",
        rm_telo="analysis_other/repeatmasker/{sample}/rm_filtered_telo.bed",
        telo="assembly/qc/{assembler}/{sample}/T2T_contigs.unphased_telomeric_repeat_windows.csv"
    output:
        gaps_all="assembly/qc/{assembler}/{sample}/karyogram.gaps.with_unassigned.png",
        gaps_clean="assembly/qc/{assembler}/{sample}/karyogram.gaps.clean.png",
        nucflag_clean="assembly/qc/{assembler}/{sample}/karyogram.nucflag.clean.png",
        gap_lengths="assembly/qc/{assembler}/{sample}/karyogram.gap_lengths.png",
        nucflag_types="assembly/qc/{assembler}/{sample}/karyogram.nucflag_types.png"
    params:
        script="workflow/scripts/plot_karyograms.R",
        output_prefix="assembly/qc/{assembler}/{sample}/karyogram"
    threads: 1
    log:
        "logs/plot_karyogram/{assembler}_{sample}_{reference}.log"
    shell:
        """
        Rscript {params.script} \
            {input.fai} \
            {input.seqinfo} \
            {input.colors} \
            {input.gap_stats} \
            {input.nucflag} \
            {input.rm_satellites} \
            {input.rm_telo} \
            {input.telo} \
            {params.output_prefix} \
            > {log} 2>&1
        """
