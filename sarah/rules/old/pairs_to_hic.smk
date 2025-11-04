# extract columns + clean strand info
rule collect_hic:
    input:
        expand("outputs/hic_files_T2T/hic/{sample}.{hap}.hic",
              sample=config["samples"],
              hap=config["hap"])

rule clean_pairs:
    input:
        "../analysis_other/porec/{sample}.{hap}/pairs/{sample}.{hap}.pairs.gz"
    output:
        "outputs/hic_files_T2T/pairs_juice/{sample}.{hap}.pairs.for_juice"
    shell:
        """
        zcat {input} | \
        grep -v '^#' | \
        awk '{{OFS="\\t"; print $1,$6,$2,$3,0,$7,$4,$5,1,$9,$10}}' | \
        awk '{{OFS="\\t"; $2=($2=="+")?0:1; $6=($6=="+")?0:1; print}}'> {output}
        """

# juicer_tools pre
rule pairs_to_hic:
    input:
        pairs="outputs/pairs_juice/{sample}.{hap}.pairs.for_juice",
        chromsize=config['ref'] + ".chrom-size.txt"
    output:
        hic = "outputs/hic_files_T2T/hic/{sample}.{hap}.hic"
    params:
        resolutions = config.get("juicer_resolutions", "1000,5000,10000,50000,100000")
    log:
        "logs/juicer_tools_pre.{sample}.{hap}.log"
    shell:
        '''
        java -Xmx50G -jar /mnt/storage2/users/ahleucs1/tools/juicertools/juicer_tools_1.19.02.jar \
        pre {input.pairs} {output.hic} {input.chromsize} -v -r {params.resolutions} > {log} 2>&1
        '''

