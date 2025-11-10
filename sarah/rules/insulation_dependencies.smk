datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_test:
    input:
        expand("outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation.tsv", 
        dataset=datasets, resolution = config['insulation_resolutions']),
        expand("outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}.png", 
        dataset=datasets, region = config['regions_fanc_names'])


rule insulation_score:
    input: 
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_nonadj_{resolution}_balanced.cool"
    output:
        insu = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation.tsv",
        bw = "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_{resolution}.insulation.tsv_bw_pseudo"   
    params:
        window = lambda wildcards: " ".join([str(mult * int(wildcards.resolution)) for mult in config['insu_window_multipliers']])
    log:
        "logs/cooltools_insulation.{dataset}_{resolution}.log"
    conda:
        "../env/cooltools.yml"
    shell:
        """
        cooltools insulation \
        --threshold Li \
        --clr-weight-name "weight" \
        --output {output.insu} \
        --bigwig \
        --verbose \
        {input.cool} {params.window} \
        >{log} 2>&1

        touch {output.bw}
        """

rule fanc_plot_imprinted_genes_new:
    input:
        insulation= "outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_25000.insulation.tsv_bw_pseudo",
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic"
    output:
        plot="outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}.png"
    params:
        insulation1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_25000.insulation.tsv.125000.bw", 
        insulation2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_25000.insulation.tsv.125000.bw",
        regions= lambda wc: wc.region.replace("_",":"),
        resolution=10000,
        name= lambda wc: wc.dataset,
        vmax=35
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.imprinted_genes.{dataset}_{region}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {output.plot} {params.regions} \
        -p mirror \
        -uvmin 0 -uvmax {params.vmax} \
        -lvmin 0 -lvmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -l hp1 hp2 -c black red --fix-chromosome {params.insulation1} {params.insulation2} \
        >{log} 2>&1
        """
