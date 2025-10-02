# Create collection of Results, Reports, and Plots for each sample
rule all_per_sample_results:
    input:
        # Global ancesty PCA plots
        expand("results/{sample}/ancestry/global_ancestry_pca.{population}.{sample}.pdf", population=["POP", "SUPERPOP"], sample=finished_samples),
        
        # Summary reports
        #"analysis_other/ancestry/reports/ancestry_summary.html"


rule plot_pca_sample:
    input:
        matrix="analysis_other/ancestry/pca/merged_cohort.eigenvec",
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