#ice correction with thresholds taken from diagnostic plot log output

rule hic_correct_matrix:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool"
        diagnostic = "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_corrected.cool"
    log:
        "logs/hic_correct.{dataset}_{exp}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),
        #filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        madscore=$(grep "mad threshold" {input.diagnostic} | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g')
        upper=$(echo -3*$madscore | bc)
        thresholds=$(echo $madscore " " $upper)

        hicCorrectMatrix correct \
        --filterThreshold $thresholds 
        --matrix {input.cool} \
        --correctionMethod {params.correction_method} \
        --outFileName {output.cool} \
        >{log} 2>&1
        """
# old:
# filter_threshold need to be set individually. but not quite feasible - what is a very conservative choice?
rule hic_correct_matrix:
    input:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}.cool",
        diagnostic = "logs/hic_diagnostic.{dataset}_{exp}.{resolution}.log"
    output:
        cool = "outputs/pairs_files_T2T/{dataset}/cooler/{dataset}_{exp}_{resolution}_corrected.cool"
    log:
        "logs/hic_correct.{dataset}_{exp}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE"),
        filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicCorrectMatrix correct \
            --matrix {input.cool} \
            --correctionMethod {params.correction_method} \
            --outFileName {output.cool} \
            --filterThreshold {params.filter_threshold} \
            >{log} 2>&1
        """
        
rule expected_cis:
    shell:
        """
        cooltools expected-cis \
        --outout {output.expected} \
        --smooth --aggregate-smoothed \
        {input.cool} \
        >{log} 2>&1
        """


#cooltools dots 
#--nproc 6 -o 'outputs/test.dots.10000.tsv' 
#--view "data/view_hg38.tsv" 
#$cool_file::resolutions/10000 
#outputs/test.expected.cis.10000.tsv

datasets = [f"{sample}.{hap}"
            for sample in config["samples"]
            for hap in config["hap"]]

rule fanc_triangle_plot:
    input:
        hic=
        insulation=
    output:
        plot=
    params:
        resolution=
    conda:
        "../env/fanc.yml"
    log:
        "logs/fanc_triangle_plot.{dataset}_{exp}_{resolution}.log"
    shell:
        """
        fancplot 
        -o {output.plot} chr18:18mb-28mb \
        -p {input.hic}@{params.resolution} \
        -m 4000000 \
        -vmin 0 -vmax 50 \
        -p line {input.insulation}
        """


