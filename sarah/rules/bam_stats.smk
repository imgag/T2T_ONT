datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule collect_qc:
    input:
        #expand("outputs/bam_stats/{dataset}/fragment_counts.stat.txt", dataset= datasets),
        #expand("outputs/bam_stats/{dataset}/fragment_counts.png", dataset= datasets),
        expand("outputs/bam_stats/{dataset}/per_read_fragment_counts_from_pairs.stat.txt", dataset= datasets),
        expand("outputs/bam_stats/{dataset}/total_fragment_counts_from_pairs.png", dataset= datasets)


rule pairtools_parse_bam:
    input:
        bam = lambda wc: f"../analysis_other/dip3d/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}/4-haplotag/{wc.dataset}.bam", 
        chromsize = config['ref'] + ".chrom-size.txt" # config['chromsize_porec']
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
        #expand_depth = config['parse2_depth']
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
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log}
        """
#--max-expansion-depth {params.expand_depth} \

rule count_reads_from_pairs:
    input:
        pairs = "outputs/pairs_files_T2T/{dataset}/pairs/{dataset}_nonadj.pairs.gz", 
    output:
        counts_per_read = "outputs/bam_stats/{dataset}/per_read_fragment_counts_from_pairs.stat.txt",
        total_counts =  "outputs/bam_stats/{dataset}/total_fragment_counts_from_pairs.stat.txt"
    log:
        "logs/count_reads_frags_from_pairs.{dataset}.log"
    shell:
        """
        zcat {input.pairs}  | grep -v "^#" | \
        awk '$8 == "UU" {{print $0}}' | \
        cut -f 1 | sort | uniq -c \
        > {output.counts_per_read}  2> {log}

        awk '{{print $1}}' {output.counts_per_read} | sort -n | uniq -c > {output.total_counts} 2>> {log}
        """

rule plot_counts:
    input:
        stats = "outputs/bam_stats/{dataset}/total_fragment_counts_from_pairs.stat.txt", 
    output:
        plot = "outputs/bam_stats/{dataset}/total_fragment_counts_from_pairs.png"
    log:
        "logs/plot_frag_hist_fragment_counts_from_pairs.{dataset}.log"
    conda: 
        "../env/r_plotting.yml"
    shell:
        """
        Rscript scripts/plot_frag_hist.R \
        --frags {input.stats} \
        --plot {output.plot} \
        >{log} 2>&1
        """





#----------------------------

