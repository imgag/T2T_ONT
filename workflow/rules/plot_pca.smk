# PCA plotting rule for Snakemake
rule plot_pca:
    input:
        matrix="data/pca/{dataset}/{file}.pca.tsv"
    output:
        plot="data/pca/{dataset}/{file}.pca_plot.pdf"
    params:
        script="workflow/scripts/plot_pca.R"
    threads: 1
    log:
        "logs/plot_pca/{dataset}_{file}.log"
    shell:
        """
        Rscript {params.script} {input.matrix} {output.plot} > {log} 2>&1
        """
