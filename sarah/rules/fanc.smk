datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_fanc:
    input:
        expand("outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}_with_gene.png", dataset=datasets, region = config['regions_fanc']),
        "data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.noexon.gtf.gz"
        #expand("outputs/pairs_files_T2T/{dataset}/plots/X_chromosome/{resolution}_{vmax}_{dataset}_X_None-None.png", dataset=datasets, resolution = [100000, 250000], vmax = [100,250]),
        #expand("outputs/pairs_files_T2T/{dataset}/plots/X_chromosome/{dataset}_{region}.png", dataset=datasets, region = ["X_100000000-140000000"])
        #expand("outputs/pairs_files_T2T/{dataset}/plots/whole_chr/{dataset}_{chrom}_None-None.png", dataset=datasets, chrom= [str(i) for i in range(1,23)] + ["X","Y"]),
        #"data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.genes.gtf.gz"


rule get_gff:
    input: 
        gtf="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gtf.gz",
        gff="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gff3.gz"
    output:
        gtf="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.noexon.gtf.gz",
        gff="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.noexon.gff3.gz"
    shell:
        """
        zcat {input.gtf}  | awk '$3 != "exon"' > {output.gtf}
        zcat {input.gff}  | awk '$3 != "exon"' > {output.gff}
        """

rule get_selected_gff:
    input: 
        gtf="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gtf.gz",
        gff="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gff3.gz"
    output:
        gtf="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.GOI.gtf.gz",
        gff="data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.GOI.gff3.gz"
    shell:
        """
        zcat {input.gtf}  | grep -E "H19|IGF2|IGF2AS|INS|INS2" > {output.gtf}
        zcat {input.gff}  | grep -E "H19|IGF2|IGF2AS|INS|INS2"  > {output.gff}
        """


rule fanc_plot_imprinted_genes_new:
    input:
        genes = "data/chm13v2.0_RefSeq_Liftoff_v5.1.sorted.gtf.gz",
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic", 
        insulation1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_25000.insulation.tsv.125000.bw", 
        insulation2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_25000.insulation.tsv.125000.bw",
        boundaries1= lambda wc: f"../analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/tad/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_25000_boundaries.gff", 
        boundaries2= lambda wc: f"../analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/tad/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_25000_boundaries.gff"
    output:
        plot="outputs/pairs_files_T2T/{dataset}/plots/imprinted_genes/{dataset}_{region}_with_gene.png"
    params:
        regions= lambda wc: wc.region.replace("_",":"),
        resolution=10000,
        name= lambda wc: wc.dataset,
        vmax=40
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_plots/fanc_triangle_plot.imprinted_genes.{dataset}_{region}.log"
    shell:
        """
        fancplot \
        --width 6 \
        -n {params.name} \
        -o {output.plot} {params.regions} \
        -p split --hide-x -vmin 0 -vmax {params.vmax} {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -l hp1 hp2 -c black red --hide-x --legend-location right --aspect-ratio 0.1 --fix-chromosome {input.insulation1} {input.insulation2} \
        -p layer --aspect-ratio 0.01 --title hp1 --hide-x --fix-chromosome {input.boundaries1} \
        -p layer --aspect-ratio 0.01 --title hp2 --hide-x --fix-chromosome {input.boundaries2} \
        -p gene -g gene_id -sq --label-field gene_id -cf black -cr red --aspect-ratio 0.3 --min-gene-size 6000 --fix-chromosome {input.genes} \
        >{log} 2>&1

        fancplot -p bar -h >>{log}
        fancplot -p layer -h >>{log}
        fancplot -p gene -h >>{log}
        """

#-p bar -l hp1 hp2 -c black red --hide-x --legend-location right --aspect-ratio 0.1 --fix-chromosome {input.boundaries1} {input.boundaries2} \
# 
# plot 1: widht 8, ar: 0.4


rule fanc_plot_X_1:
    input:
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic",
        eigs1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/compartments/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.cis.bw", 
        eigs2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/compartments/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.cis.bw"
    output:
        plot="outputs/pairs_files_T2T/{dataset}/plots/X_chromosome/{resolution}_{vmax}_{dataset}_X_None-None.png"
    params:
        regions= "X",
        resolution=lambda wc: wc.resolution,
        name= lambda wc: wc.dataset,
        vmax=lambda wc: wc.vmax
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.XCI.{dataset}_{vmax}_{resolution}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {output.plot} {params.regions} \
        -p split \
        -vmin 0 -vmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -f -l hp1 hp2 -c black red --fix-chromosome {input.eigs1} {input.eigs2} \
        >{log} 2>&1
        """

rule fanc_plot_X_2:
    input:
        hic1= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj.hic", 
        hic2= lambda wc: f"outputs/hic_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj.hic",
        insulation1= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_nonadj_40000.insulation.tsv.1000000.bw", 
        insulation2= lambda wc: f"outputs/pairs_files_T2T/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_nonadj_40000.insulation.tsv.1000000.bw"
    output:
        plot="outputs/pairs_files_T2T/{dataset}/plots/X_chromosome/{dataset}_X_100000000-140000000.png"
    params:
        regions= "X:100000000-140000000",
        resolution=40000,
        name= lambda wc: wc.dataset,
        vmax=15
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.XCI_closeup.{dataset}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {output.plot} {params.regions} \
        -p split \
        -vmin 0 -vmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -l hp1 hp2 -c black red --fix-chromosome {input.insulation1} {input.insulation2} \
        >{log} 2>&1
        """

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

