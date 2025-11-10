datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_fanc:
    input:
        #expand("outputs/pairs_files_T2T/{dataset}/plots/whole_chr/{i}_{dataset}_{chrom}_None-None.pdf", zip, dataset=datasets, i=indices_chrom,chrom=chroms),
        #expand("outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{i}_{dataset}_{region}.pdf", zip, dataset=datasets, i=indices_genes, region=imprinted_genes_filename),
        #expand("outputs/pairs_files_T2T/{dataset}/plots/insulation_chr_arms/{i}_{dataset}_{chrom}.pdf", zip, dataset=datasets, i=indices_chrom_arms, chrom = chrom_arms)
        expand("outputs/pairs_files_T2T/{dataset}/plots/whole_chr/{dataset}_{chrom}_None-None.png", dataset=datasets, chrom= [str(i) for i in range(1,23)] + ["X","Y"]),
        expand("outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}.png", dataset=datasets, config['regions_fanc'])


rule fanc_plot_whole_chr_eigs_new:
    input:
        hic= "outputs/hic_files_T2T/{dataset}/hic/{dataset}_nonadj.hic",
        eigs="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_nonadj.cis.bw"
    output:
        "outputs/pairs_files_T2T/{dataset}/plots/whole_chr/{dataset}_{chrom}_None-None.png"
    params:
        chroms= lambda wc: wc.chrom,
        resolution=250000,
        name= lambda wc: wc.dataset,
        #outdir="outputs/pairs_files_T2T/{dataset}/plots/whole_chr/",
        vmax=200
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.whole_chr.{dataset}_{chrom}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {output} {params.chroms} \
        -p square \
        -vmin 0 -vmax {params.vmax} \
        {input.hic}@{params.resolution}@KR \
        -p line -f -c black --fix-chromosome {input.eigs} \
        >{log} 2>&1
        """

rule fanc_plot_imprinted_genes_new:
    input:
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic", 
        insulation1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj_25000.insulation.tsv.125000.bw", 
        insulation2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj_25000.insulation.tsv.125000.bw",
    output:
        plot="outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}.png"
    params:
        regions= lambda wc: wc.region,
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
        -p line -l hp1 hp2 -c black red --fix-chromosome {input.insulation1} {input.insulation2} \
        >{log} 2>&1
        """

#