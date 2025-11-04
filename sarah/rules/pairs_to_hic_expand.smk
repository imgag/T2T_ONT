datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_expanded_hic:
    input:
        expand("outputs/hic_files_T2T/{dataset}/hic/{dataset}_{expansion}.hic",
              dataset=datasets,
              expansion =config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.gz",
              dataset=datasets,
              expansion =config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.stats.txt",
              dataset=datasets,
              expansion =config["expansions"]),
        expand("outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.stats.html", 
              dataset=datasets,
              expansion =config["expansions"])


rule pairtools_parse_bam:
    input:
        bam = lambda wc: f"../analysis_other/dip3d/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}/4-haplotag/{wc.dataset}.bam", 
        chromsize = config['chromsize_porec']
    output:
        pairs = "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_nonadj.pairs.gz"
    log:
        "logs/pairtools_parse2.{dataset}_nonadj.log"
    conda:
        "../env/pairtools.yml"
    params:
        assembly = "T2T-CHM13.v2",
        columns = "readID,chrom1,pos1,chrom2,pos2,strand1,strand2,mapq1,mapq2",
        orientation = "pair", 
        position = "junction",
        expand_depth = config['parse2_depth']
    shell:
        """
        pairtools --verbose parse2 \
            --chroms-path {input.chromsize} \
            --assembly {params.assembly} \
            --report-position {params.position} \
            --report-orientation {params.orientation} \
            --add-pair-index \
            --single-end \
            --expand \
            --max-expansion-depth {params.expand_depth} \
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log}
        """

# there is a -max-expansion-depth setting. should be used as there are reads ~1.6 Mb (at least in hp1+hp2 reads). 
# what cutter? what is general size of frags? I find most are <<500bp even <<250bp
# Desphande: NalIII: most reads are <21 fragments, few are <50, v few are >50. pick 20? maybe test w several

rule merge_pairs_stats:
    input:
        pairs = "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.gz"
    output:
        stats = "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.stats.txt"
    log:
        "logs/merge_pairs_stats.{dataset}_{expansion}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools stats \
            --output {output.stats} \
            {input.pairs} \
            >{log} 2>&1
        """

# Stats script from wf-pore-c pipeline
rule pairs_stats_report:
    input:
        "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.stats.txt"
    output:
        "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.stats.html"
    log:
        "logs/merge_pairs_stats_report.{dataset}_{expansion}.log"
    params:
        report_script = "scripts/create_pairs_report.py"
    conda:
        "../env/py_report.yml"
    shell:
        """
        python {params.report_script} \
            {input} {output} \
            >{log} 2>&1
        """

# extract columns + clean strand info
# readID strand1 chrom1 pos1 0 strand2 chrom2 pos2 1 mapq1 mapq2
# juicer: <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
# strand: 0 for forward, 1 for reverse
rule clean_pairs:
    input:
        "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_{expansion}.pairs.gz"
    output:
        temp("outputs/hic_files_T2T/{dataset}/pairs_juice/{dataset}_{expansion}.pairs.for_juice")
    shell:
        """
        zcat {input} | \
        grep -v '^#' | \
        awk '{{OFS="\\t"; print $1,$6,$2,$3,0,$7,$4,$5,1,$11,$12}}' | \
        awk '{{OFS="\\t"; $2=($2=="+")?0:1; $6=($6=="+")?0:1; print}}'> {output}
        """

# juicer_tools pre
rule pairs_to_hic:
    input:
        pairs = "outputs/hic_files_T2T/{dataset}/pairs_juice/{dataset}_{expansion}.pairs.for_juice",
        chromsize = config['chromsize_porec']
    output:
        hic = "outputs/hic_files_T2T/{dataset}/hic/{dataset}_{expansion}.hic"
    params:
        resolutions = config.get("juicer_resolutions", "1000,5000,10000,50000,100000")
    log:
        "logs/juicer_tools_pre/{dataset}_{expansion}.log"
    shell:
        '''
        java -Xmx50G -jar /mnt/storage2/users/ahleucs1/tools/juicertools/juicer_tools_1.19.02.jar \
        pre {input.pairs} {output.hic} {input.chromsize} -v -r {params.resolutions} > {log} 2>&1
        '''

