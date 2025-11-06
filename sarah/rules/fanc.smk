datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]


chroms = config['chroms_fanc'].split(' ')
indices_chrom = list(range(len(chroms)))

imprinted_genes_filename = [r.replace(':','_') for r in config['regions_fanc'].split(' ')]
indices_genes = list(range(len(imprinted_genes_filename)))

chrom_arms = [r.replace(':','_') for r in config['chrom_arms_fanc'].split(' ')] 
indices_chrom_arms = list(range(len(chrom_arms)))


rule collect_fanc:
    input:
        expand("outputs/pairs_files_T2T/{dataset}/plots/whole_chr/{i}_{dataset}_{chrom}_None-None.pdf", zip, dataset=datasets, i=indices_chrom,chrom=chroms),
        #expand("outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{i}_{dataset}_{region}.pdf", zip, dataset=datasets, i=indices_genes, region=imprinted_genes_filename),
        #expand("outputs/pairs_files_T2T/{dataset}/plots/insulation_chr_arms/{i}_{dataset}_{chrom}.pdf", zip, dataset=datasets, i=indices_chrom_arms, chrom = chrom_arms)


rule fanc_plot_whole_chr_eigs:
    input:
        hic= "outputs/hic_files_T2T/{dataset}/hic/{dataset}_nonadj.hic",
        eigs="outputs/pairs_files_T2T/{dataset}/compartments/{dataset}_nonadj.cis.bw"
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/whole_chr/{i}_{{dataset}}_{chrom}_None-None.pdf",
        zip,
        i=indices_chrom,
        chrom = chroms)
    params:
        chroms=config['chroms_fanc'],
        resolution=250000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/whole_chr",
        vmax=200
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.whole_chr.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.chroms} \
        -p square \
        -vmin 0 -vmax {params.vmax} \
        {input.hic}@{params.resolution}@KR \
        -p line -f -c black --fix-chromosome {input.eigs} \
        >{log} 2>&1
        """

rule fanc_plot_imprinted_genes:
    input:
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic", 
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/imprinted_genes/{i}_{{dataset}}_{region}.pdf",
        zip,
        i=indices_genes,
        region = imprinted_genes_filename)
    params:
        regions=config['regions_fanc'],
        resolution=10000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes",
        vmax=40
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.imprinted_genes.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.regions} \
        -p mirror \
        -uvmin 0 -uvmax {params.vmax} \
        -lvmin 0 -lvmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        >{log} 2>&1
        """


# -c for color, -l for log, -d max distance (pyramid)

rule fanc_plot_insulation:
    input:
        hic= "outputs/hic_files_T2T/{dataset}/hic/{dataset}_nonadj.hic",
        insulation="outputs/pairs_files_T2T/{dataset}/insulation/{dataset}_nonadj_25000.insulation.tsv.125000.bw",
    output:
        plots=expand("outputs/pairs_files_T2T/{{dataset}}/plots/insulation_chr_arms/{i}_{{dataset}}_{chrom}.pdf",
        zip,
        i=indices_chrom_arms,
        chrom = chrom_arms)
    params:
        chroms=config['chrom_arms_fanc'],
        resolution=25000,
        name= lambda wc: wc.dataset,
        outdir="outputs/pairs_files_T2T/{dataset}/plots/insulation_chr_arms",
        vmax=100
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.chrom_arms_insu.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {params.outdir} {params.chroms} \
        -p triangular \
        -vmin 0 -vmax {params.vmax} \
        {input.hic}@{params.resolution}@KR \
        -p line --fix-chromosome {input.insulation} \
        >{log} 2>&1
        """

   #-p gene --fix-chromosome {input.genes} \