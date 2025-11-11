rule all_plot:
    input:
        expand("results/{sample}/assembly_issues/{sample}.karyogram.gaps.clean.png", sample=finished_samples),
        expand("results/{sample}/ref_alignment/{sample}.synteny.png", sample=finished_samples),
        "doc/tables/for_plots/T2T_contigs.all.seqinfo.txt",
        "doc/tables/for_plots/repeatmasker_summary.all_samples.csv",
        "doc/tables/for_plots/assembly_issues_all_samples.tsv"

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

# Helper scripts to collect QC data for plotting

rule collect_T2T:
    input:
        motif_files=expand("assembly/qc/phased_verkko/{sample}/T2T_contigs.both_motif_T2T.txt", sample=finished_samples),
        alignment_files=expand("assembly/qc/phased_verkko/{sample}/T2T_contigs.both_alignment_T2T.txt", sample=finished_samples),
        gap_files=expand("assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed", sample=finished_samples)
    output:
        "doc/tables/for_plots/T2T_contigs.all.seqinfo.txt"
    params:
        samples=" ".join(finished_samples)
    log:
        "logs/pangenome/collect_T2T.log"
    shell:
        """
        python3 workflow/scripts/41_collect_T2T_contigs.py \
            --samples {params.samples} \
            --motif-files {input.motif_files} \
            --alignment-files {input.alignment_files} \
            --gap-files {input.gap_files} \
            --output {output} \
            > {log} 2>&1
        """

rule collect_repeatmasker_summaries:
    input:
        expand("analysis_other/repeatmasker/{sample}/rm_summary/{sample}_family_summary.csv", sample=finished_samples)
    output:
        "doc/tables/for_plots/repeatmasker_summary.all_samples.csv"
    params:
        samples=finished_samples
    log:
        "logs/plot/collect_repeatmasker_summaries.log"
    run:
        import pandas as pd
        
        all_data = []
        
        for i, file_path in enumerate(input):
            sample = params.samples[i]
            df = pd.read_csv(file_path)
            df['sample'] = sample
            all_data.append(df)
        
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_df = combined_df.sort_values(['sample', 'repeat_class', 'repeat_family'])
        
        # Reorder columns to put sample first
        cols = ['sample'] + [col for col in combined_df.columns if col != 'sample']
        combined_df = combined_df[cols]
        
        combined_df.to_csv(output[0], index=False)

rule collect_assembly_issues:
    input:
        flagger_files=expand("analysis_other/annotations/{sample}/flagger.lifted.bed", sample=finished_samples),
        nucflag_files=expand("analysis_other/annotations/{sample}/nucflag.lifted.bed", sample=finished_samples),
        gaps_files=expand("analysis_other/annotations/{sample}/gaps.lifted.bed", sample=finished_samples)
    output:
        "doc/tables/for_plots/assembly_issues_all_samples.tsv"
    params:
        samples=finished_samples
    log:
        "logs/plot/collect_assembly_issues.log"
    run:
        import pandas as pd
        
def parse_bed_file(file_path, issue_type, sample):
    """Parse BED file and return standardized format"""
    if not os.path.exists(file_path):
        return pd.DataFrame(columns=['chromosome', 'start', 'end', 'issue_type', 'sample'])
    
    df = pd.read_csv(file_path, sep='\t', header=None, comment='#')
    if df.empty:
        return pd.DataFrame(columns=['chromosome', 'start', 'end', 'issue_type', 'sample'])
    
    # Take first 3 columns as chr, start, end
    df = df.iloc[:, :3]
    df.columns = ['chromosome', 'start', 'end']
    
    # Convert coordinates to numeric, handling non-numeric values
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    
    # Filter out rows with invalid coordinates
    df = df.dropna(subset=['start', 'end'])
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    
    df['issue_type'] = issue_type
    df['sample'] = sample
    
    return df