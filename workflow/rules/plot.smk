rule all_plot:
    input:
        expand("results/{sample}/assembly_issues/{sample}.karyogram.gaps.clean.png", sample=finished_samples),
        expand("results/{sample}/ref_alignment/{sample}.synteny.png", sample=finished_samples),
        "doc/tables/for_plots/T2T_scaffolds.both.all.txt",
        "doc/tables/for_plots/repeatmasker_summary.all_samples.csv",
        "doc/tables/for_plots/assembly_issues_all_samples.tsv",
        "doc/tables/for_plots/cenmap_single_complete_cens.all_samples.tsv",
        "doc/tables/for_plots/cenmap_single_complete_correct_cens.all_samples.tsv",
        "doc/tables/for_plots/assembly_issues_unlifted.tsv",
        expand("doc/tables/for_plots/assembly_lengths_{dataset}.tsv", dataset=["all_samples", "hprc"]),
        expand("doc/tables/for_plots/sv_stats_details.{dataset}.tsv", dataset=["all_samples", "hprc"]),
        expand("doc/tables/for_plots/sv_stats_summary.{dataset}.tsv", dataset=["all_samples", "hprc"]),
        expand("doc/tables/for_plots/small_variant_stats_details.{dataset}.tsv", dataset=["all_samples", "hprc"]),
        expand("doc/tables/for_plots/small_variant_stats_summary.{dataset}.tsv", dataset=["all_samples", "hprc"]),
        expand("results/{sample}/ancestry/global_ancestry_pca.{population}.{sample}.pdf", population=["POP", "SUPERPOP"], sample=[s for s in finished_samples if s != "T2T_HG002"]),
        "doc/tables/for_plots/porec_stats_summary.tsv"
    output:
        "doc/tables/created_plots.done"
    shell:
        """
        touch {output}
        """

# Karyogram plotting rule for Snakemake
rule plot_karyogram:
    input:
        fai=config['ref']+".fai",
        seqinfo="assembly/qc/phased_verkko/{sample}/T2T_contigs.both.seqinfo.txt",
        colors="assembly/qc/phased_verkko/{sample}/colors.tsv",
        gap_stats="assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed",
        nucflag="analysis_other/nucflag/{sample}/nucflag_misasm.bed",
        rm_satellites="analysis_other/repeatmasker/{sample}/rm_summary/{sample}_filtered_satellites.bed",
        telo="assembly/qc/phased_verkko/{sample}/T2T_contigs.both_motif.csv"
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

rule plot_pca_sample:
    input:
        matrix="analysis_other/ancestry/pca/merge_ref_samples.eigenvec",
        metadata="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam",
    output:
        plot="results/{sample}/ancestry/global_ancestry_pca.{population}.{sample}.pdf"
    params:
        script="workflow/scripts/30_plot_pca.R"
    threads: 1
    log:
        "logs/per_sample_results/{sample}/pca_plot_{population}.log"
    conda:
        "../env/r_ancestry.yml"
    shell:
        """
        Rscript \
            {params.script} \
            {input.matrix} \
            {output.plot} \
            {wildcards.sample} \
            {input.metadata} \
            {wildcards.population} > {log} 2>&1
        """


# Helper scripts to collect QC data for plotting

rule collect_T2T:
    input:
        seginfo_files=expand("assembly/qc/phased_verkko/{sample}/T2T_scaffolds.both.txt", sample=finished_samples),
    output:
        "doc/tables/for_plots/T2T_scaffolds.both.all.txt"
    log:
        "logs/plot/collect_T2T.log"
    shell:
        """
        head -n 1 $(echo {input.seginfo_files} | cut -d' ' -f1) > {output} 2> {log}
        tail -n +2 -q {input.seginfo_files} >> {output} 2>> {log}
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

rule collect_unlifted_assembly_issues:
    input:
        flagger=expand("analysis_other/flagger/{sample}/final_flagger_prediction.bed", sample=finished_samples),
        nucflag=expand("analysis_other/nucflag/{sample}/nucflag_misasm.bed", sample=finished_samples),
        gaps=expand("assembly/qc/phased_verkko/{sample}/gap_stats.both.n_regions.bed", sample=finished_samples)
    output:
        "doc/tables/for_plots/assembly_issues_unlifted.tsv"
    params:
        samples=finished_samples
    conda:
        "../env/py_report.yml"
    log:
        "logs/plot/collect_unlifted_assembly_issues.log"
    shell:
        """
        python3 workflow/scripts/55_merge_raw_issue_beds.py \
            --flagger {input.flagger} \
            --nucflag {input.nucflag} \
            --gaps {input.gaps} \
            --samples {params.samples} \
            --output {output} \
            > {log} 2>&1
        """

rule collect_assembly_issues:
    input:
        flagger_files=expand("analysis_other/annotations/{sample}/flagger.lifted.bed", sample=finished_samples),
        nucflag_files=expand("analysis_other/annotations/{sample}/nucflag.lifted.bed", sample=finished_samples),
        gaps_files=expand("analysis_other/annotations/{sample}/gaps.lifted.bed", sample=finished_samples)
    output:
        "doc/tables/for_plots/assembly_issues_all_samples.tsv"
    params:
        samples=finished_samples,
        n_samples=len(finished_samples)
    log:
        "logs/plot/collect_assembly_issues.log"
    shell:
        """
        echo -e "sample\tsource\tchr\tstart\tend\tfeature_type" > {output} 2> {log}
        
        for i in $(seq 0 $(({params.n_samples}-1))); do
            sample=$(echo "{params.samples}" | cut -d' ' -f$((i+1)))
            flagger=$(echo "{input.flagger_files}" | cut -d' ' -f$((i+1)))
            nucflag=$(echo "{input.nucflag_files}" | cut -d' ' -f$((i+1)))
            gaps=$(echo "{input.gaps_files}" | cut -d' ' -f$((i+1)))
            
            awk -v s="$sample" -v t="Flagger" 'BEGIN{{{{OFS="\\t"}}}} !/^#/ && NF>=3 {{{{print s,t,$1,$2,$3,(NF>=4? $4:"unknown")}}}}' "$flagger" >> {output} 2>> {log}
            awk -v s="$sample" -v t="Nucflag" 'BEGIN{{{{OFS="\\t"}}}} !/^#/ && NF>=3 {{{{print s,t,$1,$2,$3,(NF>=4?$4:"unknown")}}}}' "$nucflag" >> {output} 2>> {log}
            awk -v s="$sample" -v t="GapAnalysis" 'BEGIN{{{{OFS="\\t"}}}} !/^#/ && NF>=3 {{{{print s,t,$1,$2,$3,(NF>=4?$4:"unknown")}}}}' "$gaps" >> {output} 2>> {log}
        done
        """

rule collect_cenmap_single:
    input:
        bed_c = expand("analysis_other/cenmap_single/{asm}/results/final/bed/{asm}_complete_cens.bed", asm = finished_samples),
        bed_cc = expand("analysis_other/cenmap_single/{asm}/results/final/bed/{asm}_complete_correct_cens.bed", asm = finished_samples),
    output:
        tsv_c = "doc/tables/for_plots/cenmap_single_complete_cens.all_samples.tsv",
        tsv_cc = "doc/tables/for_plots/cenmap_single_complete_correct_cens.all_samples.tsv"
    params:
            script = "workflow/scripts/47_merge_cenmap_beds.py"
    log:
        "logs/plot/collect_cenmap_single.log"
    shell:
        """
        python {params.script} -i {input.bed_c} -o {output.tsv_c} 2> {log}
        python {params.script} -i {input.bed_cc} -o {output.tsv_cc} 2>> {log}
        """

hprc_samples, = glob_wildcards("data/ref/hprc/{sample}/{sample}.haplotype1.fasta")

def get_assembly_inputs(wildcards):
    """Get appropriate input files based on dataset wildcard"""
    if wildcards.dataset == "all_samples":
        return {
            "fai_hap1": expand("assembly/output/verkko/{sample}/assembly.haplotype1.fasta.fai", sample=finished_samples),
            "fai_hap2": expand("assembly/output/verkko/{sample}/assembly.haplotype2.fasta.fai", sample=finished_samples)
        }
    elif wildcards.dataset == "hprc":
        return {
            "fai_hap1": expand("data/ref/hprc/{sample}/{sample}.haplotype1.fasta.fai", sample=hprc_samples),
            "fai_hap2": expand("data/ref/hprc/{sample}/{sample}.haplotype2.fasta.fai", sample=hprc_samples)
        }

def get_assembly_length_samples(wildcards):
    """Get appropriate sample list based on dataset wildcard"""
    if wildcards.dataset == "all_samples":
        return finished_samples
    elif wildcards.dataset == "hprc":
        return hprc_samples


def get_sv_stats_inputs(wildcards):
    """Get SV stats input files based on dataset wildcard."""
    if wildcards.dataset == "all_samples":
        return {
            "details": expand("assembly/qc/phased_verkko/{sample}/sv_stats_details.tsv", sample=finished_samples),
            "summaries": expand("assembly/qc/phased_verkko/{sample}/sv_stats_summary.tsv", sample=finished_samples)
        }
    elif wildcards.dataset == "hprc":
        return {
            "details": expand("assembly/qc/phased_hprc/{sample}/sv_stats_details.tsv", sample=hprc_samples),
            "summaries": expand("assembly/qc/phased_hprc/{sample}/sv_stats_summary.tsv", sample=hprc_samples)
        }


def get_sv_stats_samples(wildcards):
    """Get sample list for SV stats based on dataset wildcard."""
    if wildcards.dataset == "all_samples":
        return finished_samples
    elif wildcards.dataset == "hprc":
        return hprc_samples


def get_small_variant_stats_inputs(wildcards):
    """Get small variant stats input files based on dataset wildcard."""
    if wildcards.dataset == "all_samples":
        return {
            "details": expand("assembly/qc/phased_verkko/{sample}/small_variant_stats_details.tsv", sample=finished_samples),
            "summaries": expand("assembly/qc/phased_verkko/{sample}/small_variant_stats_summary.tsv", sample=finished_samples)
        }
    elif wildcards.dataset == "hprc":
        return {
            "details": expand("assembly/qc/phased_hprc/{sample}/small_variant_stats_details.tsv", sample=hprc_samples),
            "summaries": expand("assembly/qc/phased_hprc/{sample}/small_variant_stats_summary.tsv", sample=hprc_samples)
        }


def get_small_variant_stats_samples(wildcards):
    """Get sample list for small variant stats based on dataset wildcard."""
    if wildcards.dataset == "all_samples":
        return finished_samples
    elif wildcards.dataset == "hprc":
        return hprc_samples

rule collect_assembly_lengths:
    input:
        unpack(get_assembly_inputs)
    output:
        "doc/tables/for_plots/assembly_lengths_{dataset}.tsv"
    params:
        samples=get_assembly_length_samples
    log:
        "logs/plot/collect_assembly_lengths_{dataset}.log"
    run:
        import pandas as pd
        
        all_data = []
        
        # Process haplotype1
        for i, fai_file in enumerate(input.fai_hap1):
            sample = params.samples[i]
            df = pd.read_csv(fai_file, sep='\t', header=None, 
                           names=['contig', 'length', 'offset', 'linebases', 'linewidth'])
            df['sample'] = sample
            df['haplotype'] = 'haplotype1'
            all_data.append(df[['sample', 'haplotype', 'contig', 'length']])
        
        # Process haplotype2
        for i, fai_file in enumerate(input.fai_hap2):
            sample = params.samples[i]
            df = pd.read_csv(fai_file, sep='\t', header=None,
                           names=['contig', 'length', 'offset', 'linebases', 'linewidth'])
            df['sample'] = sample
            df['haplotype'] = 'haplotype2'
            all_data.append(df[['sample', 'haplotype', 'contig', 'length']])
        
        # Combine all data
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_df = combined_df.sort_values(['sample', 'haplotype', 'length'], ascending=[True, True, False])
        
        # Save to output
        combined_df.to_csv(output[0], sep='\t', index=False)


rule collect_sv_stats_details:
    """Collect detailed SV statistics from samples (dataset-specific)."""
    input:
        unpack(get_sv_stats_inputs)
    output:
        "doc/tables/for_plots/sv_stats_details.{dataset}.tsv"
    params:
        samples=get_sv_stats_samples
    log:
        "logs/plot/collect_sv_stats_details_{dataset}.log"
    run:
        import pandas as pd
        
        all_data = []
        if input.details:
            for f in input.details:
                df = pd.read_csv(f, sep='\t')
                all_data.append(df)
        
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            combined = combined.sort_values(['sample', 'chrom', 'pos'])
            combined.to_csv(output[0], sep='\t', index=False)
        else:
            # Create empty file with header
            pd.DataFrame(columns=['sample', 'chrom', 'pos', 'end', 'sv_id', 'sv_type', 
                                  'sv_len', 'size_class', 'genotype', 'is_phased', 
                                  'zygosity', 'haplotype', 'filter']).to_csv(output[0], sep='\t', index=False)


rule collect_sv_stats_summary:
    """Collect summary SV statistics from samples (dataset-specific)."""
    input:
        unpack(get_sv_stats_inputs)
    output:
        "doc/tables/for_plots/sv_stats_summary.{dataset}.tsv"
    params:
        samples=get_sv_stats_samples
    log:
        "logs/plot/collect_sv_stats_summary_{dataset}.log"
    run:
        import pandas as pd
        
        all_data = []
        if input.summaries:
            for f in input.summaries:
                df = pd.read_csv(f, sep='\t')
                all_data.append(df)
        
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            combined = combined.sort_values(['sample', 'category', 'sv_type', 'subcategory'])
            combined.to_csv(output[0], sep='\t', index=False)
        else:
            pd.DataFrame(columns=['sample', 'category', 'subcategory', 'sv_type', 'count', 
                                  'total_length', 'mean_length']).to_csv(output[0], sep='\t', index=False)


rule collect_small_variant_stats_details:
    """Collect detailed small variant statistics from samples (dataset-specific)."""
    input:
        unpack(get_small_variant_stats_inputs)
    output:
        "doc/tables/for_plots/small_variant_stats_details.{dataset}.tsv"
    params:
        samples=get_small_variant_stats_samples
    log:
        "logs/plot/collect_small_variant_stats_details_{dataset}.log"
    run:
        import pandas as pd
        
        all_data = []
        if input.details:
            for f in input.details:
                df = pd.read_csv(f, sep='\t')
                all_data.append(df)
        
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            combined = combined.sort_values(['sample', 'chrom', 'pos'])
            combined.to_csv(output[0], sep='\t', index=False)
        else:
            pd.DataFrame(columns=['sample', 'chrom', 'pos', 'var_id', 'ref', 'alt', 
                                  'var_type', 'var_len', 'genotype', 'is_phased', 
                                  'zygosity', 'haplotype', 'filter', 'ti_tv', 
                                  'is_biallelic', 'allele_freq']).to_csv(output[0], sep='\t', index=False)


rule collect_porec_stats:
    """Aggregate pairtools pairs stats across all samples (unphased + hp1/hp2).

    Three files per finished sample are collected:
      - analysis_other/porec/{sample}/          (unphased, wf-pore-c)
      - analysis_other/porec/{sample}.hp1/      (haplotype 1, dip3d)
      - analysis_other/porec/{sample}.hp2/      (haplotype 2, dip3d)
    """
    input:
        stats=expand(
            "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.txt",
            dataset=[f"{s}{hp}" for s in finished_samples for hp in ["", ".hp1", ".hp2"]]
        )
    output:
        "doc/tables/for_plots/porec_stats_summary.tsv"
    params:
        script="workflow/scripts/58_aggregate_porec_stats.py"
    conda:
        "../env/py_report.yml"
    log:
        "logs/plot/collect_porec_stats.log"
    shell:
        """
        python3 {params.script} \\
            --stats {input.stats} \\
            --output {output} \\
            > {log} 2>&1
        """


rule collect_small_variant_stats_summary:
    """Collect summary small variant statistics from samples (dataset-specific)."""
    input:
        unpack(get_small_variant_stats_inputs)
    output:
        "doc/tables/for_plots/small_variant_stats_summary.{dataset}.tsv"
    params:
        samples=get_small_variant_stats_samples
    log:
        "logs/plot/collect_small_variant_stats_summary_{dataset}.log"
    run:
        import pandas as pd
        
        all_data = []
        if input.summaries:
            for f in input.summaries:
                df = pd.read_csv(f, sep='\t')
                all_data.append(df)
        
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            combined = combined.sort_values(['sample', 'category', 'var_type', 'subcategory'])
            combined.to_csv(output[0], sep='\t', index=False)
        else:
            pd.DataFrame(columns=['sample', 'category', 'subcategory', 'var_type',
                                  'count', 'percentage']).to_csv(output[0], sep='\t', index=False)


# ── Assembly issue classification (all samples) ───────────────────────────────
#
# Three-way split of NucFlag / Flagger issue regions by repeat context:
#   satellite      RepeatMasker Satellite/* (alpha-sat, HSAT, acro, subtelo)
#   other_repeat   All remaining annotated repeats (LINE, SINE, LTR, …)
#   non_repetitive No RepeatMasker annotation
#
# All coordinates are in Verkko assembly space; RepeatMasker .out is used
# directly — no liftover needed.

_ISSUE_TYPES = {
    "nucflag": ["COLLAPSE", "COLLAPSE_OTHER", "COLLAPSE_VAR", "HET", "MISJOIN"],
    "flagger": ["Col", "Dup", "Err", "Hap"],
}
_ALL_TOOL_ISSUES = [(t, i) for t, issues in _ISSUE_TYPES.items() for i in issues]
_ALL_ISSUES      = [i for _, i in _ALL_TOOL_ISSUES]


rule all_assembly_issue_classification:
    input:
        expand(
            "analysis_other/assembly_issue_classification/{sample}/issue_summary.tsv",
            sample=finished_samples,
        ),


rule rm_satellite_bed:
    """Merge all Satellite/* entries from RepeatMasker into a single BED."""
    input:
        "analysis_other/repeatmasker/{sample}/assembly.fasta.out",
    output:
        "analysis_other/assembly_issue_classification/{sample}/repeat_annot/satellite.bed",
    log:
        "logs/assembly_issue_classification/{sample}/rm_satellite_bed.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        awk 'NR>2 && $11~/^Satellite/ {{print $5"\t"($6-1)"\t"$7}}' {input} \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > {output} 2>{log}
        """


rule rm_all_repeat_bed:
    """Merge all annotated repeat entries from RepeatMasker into a single BED."""
    input:
        "analysis_other/repeatmasker/{sample}/assembly.fasta.out",
    output:
        "analysis_other/assembly_issue_classification/{sample}/repeat_annot/all_repeat.bed",
    log:
        "logs/assembly_issue_classification/{sample}/rm_all_repeat_bed.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        awk 'NR>2 && $11!="" {{print $5"\t"($6-1)"\t"$7}}' {input} \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            > {output} 2>{log}
        """


rule nucflag_issue_bed:
    """Extract one NucFlag issue type into a 3-column BED."""
    input:
        "analysis_other/nucflag/{sample}/nucflag_misasm.bed",
    output:
        "analysis_other/assembly_issue_classification/{sample}/issues/nucflag/{issue}.bed",
    wildcard_constraints:
        issue = "|".join(_ISSUE_TYPES["nucflag"]),
    log:
        "logs/assembly_issue_classification/{sample}/nucflag_{issue}.log",
    shell:
        """
        mkdir -p $(dirname {output})
        awk '$4=="{wildcards.issue}"' {input} \
            | cut -f1-3 | sort -k1,1 -k2,2n \
            > {output} 2>{log}
        """


rule flagger_issue_bed:
    """Extract one Flagger issue type into a 3-column BED."""
    input:
        "analysis_other/flagger/{sample}/final_flagger_prediction.bed",
    output:
        "analysis_other/assembly_issue_classification/{sample}/issues/flagger/{issue}.bed",
    wildcard_constraints:
        issue = "|".join(_ISSUE_TYPES["flagger"]),
    log:
        "logs/assembly_issue_classification/{sample}/flagger_{issue}.log",
    shell:
        """
        grep -v '^track' {input} \
            | awk '$4=="{wildcards.issue}"' \
            | cut -f1-3 | sort -k1,1 -k2,2n \
            > {output} 2>{log}
        """


rule issue_intersect_satellite:
    """Issue regions overlapping satellite repeats."""
    input:
        issues    = "analysis_other/assembly_issue_classification/{sample}/issues/{tool}/{issue}.bed",
        satellite = "analysis_other/assembly_issue_classification/{sample}/repeat_annot/satellite.bed",
    output:
        "analysis_other/assembly_issue_classification/{sample}/intersect/{tool}/{issue}_satellite.bed",
    wildcard_constraints:
        tool  = "nucflag|flagger",
        issue = "|".join(_ALL_ISSUES),
    log:
        "logs/assembly_issue_classification/{sample}/{tool}_{issue}_satellite.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        mkdir -p $(dirname {output})
        bedtools intersect -a {input.issues} -b {input.satellite} \
            | sort -k1,1 -k2,2n | bedtools merge -i - \
            > {output} 2>{log}
        """


rule issue_intersect_other_repeat:
    """Issue regions in non-satellite annotated repeats (all_repeat minus satellite)."""
    input:
        issues     = "analysis_other/assembly_issue_classification/{sample}/issues/{tool}/{issue}.bed",
        all_repeat = "analysis_other/assembly_issue_classification/{sample}/repeat_annot/all_repeat.bed",
        satellite  = "analysis_other/assembly_issue_classification/{sample}/repeat_annot/satellite.bed",
    output:
        "analysis_other/assembly_issue_classification/{sample}/intersect/{tool}/{issue}_other_repeat.bed",
    wildcard_constraints:
        tool  = "nucflag|flagger",
        issue = "|".join(_ALL_ISSUES),
    log:
        "logs/assembly_issue_classification/{sample}/{tool}_{issue}_other_repeat.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        bedtools intersect -a {input.issues} -b {input.all_repeat} \
            | bedtools subtract -a - -b {input.satellite} \
            | sort -k1,1 -k2,2n | bedtools merge -i - \
            > {output} 2>{log}
        """


rule issue_intersect_non_repetitive:
    """Issue regions with no RepeatMasker annotation."""
    input:
        issues     = "analysis_other/assembly_issue_classification/{sample}/issues/{tool}/{issue}.bed",
        all_repeat = "analysis_other/assembly_issue_classification/{sample}/repeat_annot/all_repeat.bed",
    output:
        "analysis_other/assembly_issue_classification/{sample}/intersect/{tool}/{issue}_non_repetitive.bed",
    wildcard_constraints:
        tool  = "nucflag|flagger",
        issue = "|".join(_ALL_ISSUES),
    log:
        "logs/assembly_issue_classification/{sample}/{tool}_{issue}_non_repetitive.log",
    conda:
        "../env/gqc.yml"
    shell:
        """
        bedtools subtract -a {input.issues} -b {input.all_repeat} \
            > {output} 2>{log}
        """


rule assembly_issue_classification:
    """
    Aggregate all issue × repeat-class intersections into:
      issue_summary.tsv  — aggregate bp / region counts per tool × issue × repeat class
      issue_regions.tsv  — per-region rows for size histograms
    """
    input:
        intersect = lambda wc: [
            f"analysis_other/assembly_issue_classification/{wc.sample}"
            f"/intersect/{tool}/{issue}_{rclass}.bed"
            for tool, issue in _ALL_TOOL_ISSUES
            for rclass in ["satellite", "other_repeat", "non_repetitive"]
        ],
        issues = lambda wc: [
            f"analysis_other/assembly_issue_classification/{wc.sample}"
            f"/issues/{tool}/{issue}.bed"
            for tool, issue in _ALL_TOOL_ISSUES
        ],
    output:
        summary = "analysis_other/assembly_issue_classification/{sample}/issue_summary.tsv",
        regions = "analysis_other/assembly_issue_classification/{sample}/issue_regions.tsv",
    log:
        "logs/assembly_issue_classification/{sample}/assembly_issue_classification.log",
    run:
        import pandas as pd

        log_fh  = open(log[0], "w")
        sample  = wildcards.sample
        rclasses = ["satellite", "other_repeat", "non_repetitive"]

        def read_bed(path):
            rows = []
            with open(path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    f = line.split("\t")
                    rows.append((f[0], int(f[1]), int(f[2])))
            return rows

        summary_rows = []
        region_rows  = []

        for tool, issue in _ALL_TOOL_ISSUES:
            all_regions = read_bed(
                f"analysis_other/assembly_issue_classification/{sample}"
                f"/issues/{tool}/{issue}.bed"
            )
            total_bp = sum(e - s for _, s, e in all_regions)
            print(f"{tool}/{issue}: {len(all_regions)} regions, {total_bp:,} bp", file=log_fh)

            for rclass in rclasses:
                path = (
                    f"analysis_other/assembly_issue_classification/{sample}"
                    f"/intersect/{tool}/{issue}_{rclass}.bed"
                )
                regions = read_bed(path)
                bp  = sum(e - s for _, s, e in regions)
                pct = round(100 * bp / total_bp, 2) if total_bp else 0
                summary_rows.append(dict(
                    sample       = sample,
                    tool         = tool,
                    issue_type   = issue,
                    repeat_class = rclass,
                    n_regions    = len(regions),
                    total_bp     = bp,
                    pct_of_total = pct,
                ))
                for chrom, start, end in regions:
                    region_rows.append(dict(
                        sample       = sample,
                        tool         = tool,
                        issue_type   = issue,
                        repeat_class = rclass,
                        chrom        = chrom,
                        start        = start,
                        end          = end,
                        length_bp    = end - start,
                    ))
                print(f"  {rclass}: {len(regions)} regions, {bp:,} bp ({pct}%)", file=log_fh)

        pd.DataFrame(summary_rows).to_csv(output.summary, sep="\t", index=False)
        pd.DataFrame(region_rows).to_csv(output.regions,  sep="\t", index=False)
        print("Done.", file=log_fh)
        log_fh.close()


rule plot_Y_karyogram:
    """Y chromosome length + error-track overview across all finished samples."""
    input:
        samples_yml = "data/finished_samples.yml",
        colors       = expand("assembly/qc/phased_verkko/{sample}/colors.tsv",
                              sample=finished_samples),
        paf          = expand("assembly/qc/phased_verkko/{sample}/both.mapped_T2T.paf",
                              sample=finished_samples),
    output:
        pdf      = "doc/figures/fig_Y_karyogram.Y_karyogram.pdf",
        png      = "doc/figures/fig_Y_karyogram.Y_karyogram.png",
        frag_pdf = "doc/figures/fig_Y_karyogram.Y_karyogram_fragmented.pdf",
        frag_png = "doc/figures/fig_Y_karyogram.Y_karyogram_fragmented.png",
    params:
        script        = "workflow/scripts/59_plot_Y_karyogram.R",
        output_prefix = "doc/figures/fig_Y_karyogram",
    conda:
        "../env/R.yml"
    log:
        "logs/plot/Y_karyogram.log"
    shell:
        """
        Rscript {params.script} \
            {input.samples_yml} \
            {params.output_prefix} \
            > {log} 2>&1
        """

