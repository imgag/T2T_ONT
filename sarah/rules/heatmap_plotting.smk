datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]


#hicExplorer Plotting Functions
#todo: hicPlotTADs: pyramid plot with tracks. but; weird .ini file necessary; how to implement with snakemake? can i dynamically write the ini?

rule collect_plot:
    input:
        expand("outputs/pairs_files_T2T/{dataset}/plots/{dataset}_{exp}_{resolution}_with_eigs.png",
        dataset=datasets,
        exp=config["expansions"],
        resolution=[100000]),
        expand("outputs/pairs_files_T2T/{dataset}/plots/{dataset}_{exp}_{resolution}_with_insulation.png",
        dataset=datasets,
        exp=config["expansions"],
        resolution=[25000])

rule plot_matrix_eigs:
    input: 
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_balanced.cool",
        eigs = "outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_{exp}.cis.bw"
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/plots/{dataset}_{exp}_{resolution}_with_eigs.png"
    params:
        region = "chr1:0-245000000"
    conda:
        "../env/hicexplorer.yml"
    log:
        "logs/hicPlotMatrix_eigs.{dataset}_{exp}_{resolution}.log"
    shell:
        """
        hicPlotMatrix \
        --matrix {input.cool} \
        --outFileName {output.plot} \
        --region {params.region} \
        --vMax 0.01 \
        --log \
        --bigwig {input.eigs} \
        >{log} 2>&1
        """

rule plot_matrix_insulation:
    input: 
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_balanced.cool",
        insu = expand("outputs/pairs_files_T2T/{{dataset}}/insulation/{{dataset}}_{{exp}}_{{resolution}}.insulation.tsv.{window}.bw",
        #window=[120000,200000,400000,1000000] # for 40 kb
        window=[75000,125000,250000,625000] # for 25 kb
        )
    output:
        plot = "outputs/pairs_files_T2T/{dataset}/plots/{dataset}_{exp}_{resolution}_with_insulation.png"
    params:
        region = "chr1:0-60000000"
    conda:
        "../env/hicexplorer.yml"
    log:
        "logs/hicPlotMatrix_insulation.{dataset}_{exp}_{resolution}.log"
    shell:
        """
        hicPlotMatrix \
        --matrix {input.cool} \
        --outFileName {output.plot} \
        --region {params.region} \
        --vMax 0.01 \
        --log \
        --bigwig {input.insu} \
        >{log} 2>&1
        """


