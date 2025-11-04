
rule collect_pairs_files:
    input:
        expand("outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.gz",
              sample=config["samples"],
              hap=config["hap"],
              expansion = ["nonadj","adj"]),
        expand("outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.stats.txt",
              sample=config["samples"],
              hap=config["hap"],
              expansion = ["nonadj","adj"]),
        expand("outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.stats.html", 
              sample=config["samples"], 
              hap=config["hap"],
              expansion = ["nonadj","adj"])

rule pairtools_parse_bam_expand:
    input:
        bam = "../analysis_other/dip3d/{sample}/4-haplotag/{sample}.{hap}.bam", 
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        pairs = "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_nonadj.pairs.gz"
    log:
        "logs/pairtools_parse2.{sample}.{hap}_nonadj.log"
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

rule pairtools_parse_bam_no_expand:
    input:
        bam = "../analysis_other/dip3d/{sample}/4-haplotag/{sample}.{hap}.bam", 
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        pairs = "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_adj.pairs.gz"
    log:
        "logs/pairtools_parse2_no_expand.{sample}.{hap}_adj.log"
    conda:
        "../env/pairtools.yml"
    params:
        assembly = "T2T-CHM13.v2",
        columns = "readID,chrom1,pos1,chrom2,pos2,strand1,strand2,mapq1,mapq2",
        orientation = "pair", 
        position = "junction"
    shell:
        """
        pairtools --verbose parse2 \
            --chroms-path {input.chromsize} \
            --assembly {params.assembly} \
            --report-position {params.position} \
            --report-orientation {params.orientation} \
            --add-pair-index \
            --single-end \
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log}
        """

rule merge_pairs_stats:
    input:
        pairs = "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.gz"
    output:
        stats = "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.stats.txt"
    log:
        "logs/merge_pairs_stats.{sample}.{hap}_{expansion}.log"
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
        "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.stats.txt"
    output:
        "outputs/pairs_files_T2T/{sample}.{hap}/pairs/{sample}.{hap}_{expansion}.pairs.stats.html"
    log:
        "logs/merge_pairs_stats_report.{sample}.{hap}_{expansion}.log"
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
