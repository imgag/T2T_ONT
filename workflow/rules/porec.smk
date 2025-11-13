rule all_porec:
    input:
        expand("analysis_other/porec/{sample}.done", sample = "T2T01")

rule collect_porec:
    input:
        expand("analysis_other/porec/{{asm}}{hp}/pairs/{{asm}}{hp}.pairs.gz", hp = ["", ".hp1", ".hp2"]),
        #expand("analysis_other/porec/{{asm}}{hp}/pairs/{{asm}}{hp}_adj.pairs.gz", hp = ["", ".hp1", ".hp2"]), # new
        expand("analysis_other/porec/{{asm}}{hp}/pairs/{{asm}}{hp}.pairs.stats.html", hp = ["", ".hp1", ".hp2"]),
        expand("analysis_other/porec/{{asm}}{hp}/cooler/{{asm}}{hp}.mcool", hp = ["", ".hp1", ".hp2"]),
        expand("analysis_other/porec/{{asm}}{hp}/cooler/{{asm}}{hp}_{res}_corrected.cool", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/qc/{{asm}}{hp}_{res}_diagnostic.png", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/qc/plot_vs_counts_{res}.png", hp = ["", ".hp1", ".hp2"], res=config['porec_resolutions']),
        expand("analysis_other/porec/{{asm}}{hp}/tad/{{asm}}{hp}_{res}_domains.bed", hp = ["", ".hp1", ".hp2"], res=config.get('tad_resolutions', ['25000'])),
        expand("analysis_other/porec/{{asm}}{hp}/loops/{{asm}}{hp}_{res}_loops.bedpe", hp = ["", ".hp1", ".hp2"], res=config.get('loop_resolutions', ['10000'])),
        expand("analysis_other/porec/{{asm}}{hp}/hic/{{asm}}{hp}.hic", hp = ["", ".hp1", ".hp2"]), # new
        expand("analysis_other/porec/{{asm}}{hp}/compartments/{{asm}}{hp}.cis{file}", hp = ["", ".hp1", ".hp2"] ,file=[".vecs.tsv",".lam.txt",".bw"]),  # new
        expand("analysis_other/porec/{{asm}}{hp}/insulation/{{asm}}{hp}_{res}.insulation.tsv", hp = ["", ".hp1", ".hp2"], res = config['insulation_resolutions']), # new
        #expand("analysis_other/porec/{{asm}}{hp}/plots/whole_chr/{{asm}}{hp}_{chrom}_None-None.png", hp = ["", ".hp1", ".hp2"], chrom= [str(i) for i in range(1,23)] + ["X", "Y"] ), # new
        #expand("analysis_other/porec/{{asm}}{hp}/plots/imprinted_genes/{{asm}}{hp}_{region}.png", hp = [".hp1", ".hp2"], region = config['regions_fanc']) # new
    output:
        "analysis_other/porec/{asm}.done"
    shell:
        """
        touch {output}
        """

def get_all_porec_runs(wc):
    folders = datasets[wc.dataset].get("POREC", [])
    # If folders is not a list, make it one
    if not isinstance(folders, list):
        folders = [folders] if folders else []
    
    run_ids = []
    for folder_path in folders:
        # Extract the basename (run ID) from the folder path
        run_id = os.path.basename(folder_path)
        # Look up the run_id in unique_datasets["porec"] to verify it exists
        if run_id in unique_datasets["porec"]:
            run_ids.append(run_id)

    return run_ids

# Create pairs file for Haplotype1 and Haplotype2 datasets
# 1. adj contacts (_adj, not used downstream)
rule pairtools_parse_bam_no_expand:
    input:
        bam = lambda wc: f"analysis_other/dip3d/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}/4-haplotag/{wc.dataset}.tagged.  bam", 
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}_adj.pairs.gz"
    log:
        "logs/porec/pairtools_parse2_no_expand.{dataset}_adj.log"
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
        pairtools parse2 \
            --chroms-path {input.chromsize} \
            --assembly {params.assembly} \
            --report-position {params.position} \
            --report-orientation {params.orientation} \
            --single-end \
            --readid-transform 'readID.split(":")[0]' \
            --drop-seq \
            --drop-sam \
            --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,mismatches \
            --output {output.pairs} \
            {input.bam} \
            >{log}
        """

# 2. adj+nonadj contacts (no suffix, used downstream)
rule pairtools_parse_bam:
    input:
        bam = lambda wc: f"analysis_other/dip3d/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}/4-haplotag/{wc.dataset}.bam", 
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    log:
        "logs/porec/pairtools_parse2.{dataset}.log"
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
        pairtools parse2 \
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

rule pairtools_parse_bam_wfporec:
    input:
        chromsize = config['ref'] + ".chrom-size.txt",
        bam = "analysis_other/wf-pore-c/{run}/bams/{run}.ns.bam"
    output:
        pairs = "analysis_other/wf-pore-c/{run}/pairs/{run}.nonadj.pairs.gz"
    log:
        "logs/porec/pairtools_parse2_nonadj.{run}.log"
    conda:
        "../env/pairtools.yml"  
    params:
        assembly = "T2T-CHM13.v2",
        columns = "readID,chrom1,pos1,chrom2,pos2,strand1,strand2,mapq1,mapq2",
        orientation = "pair", 
        position = "junction",
        expand_depth = config['parse2_depth']
    threads:
        12
    shell:
        """
        pairtools parse2 \
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
            --nproc-in 3 \
            --nproc-out {threads} \
            {input.bam} \
            >{log}
        """

# Create pairs file for all flowcells from sample
rule merge_pairs:
    input:
        pairs = lambda wc: expand("analysis_other/wf-pore-c/{run}/pairs/{run}.nonadj.pairs.gz", run = get_all_porec_runs(wc))
    output:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    wildcard_constraints:
        dataset=r"\w+"
    log:
        "logs/porec/merge_pairs.{dataset}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools merge \
            --output {output.pairs} \
            {input.pairs} \
            >{log} 2>&1
        """

rule merge_pairs_stats:
    input:
        stats = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    output:
        stats = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.txt"
    log:
        "logs/porec/merge_pairs_stats.{dataset}.log"
    conda:
        "../env/pairtools.yml"
    shell:
        """
        pairtools stats \
            --output {output.stats} \
            {input.stats} \
            >{log} 2>&1
        """

# Stats script from wf-pore-c pipeline
rule pairs_stats_report:
    input:
        "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.txt"
    output:
        "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.stats.html"
    log:
        "logs/porec/merge_pairs_stats_report.{dataset}.log"
    params:
        report_script = "workflow/scripts/create_pairs_report.py"
    conda:
        "../env/py_report.yml"
    shell:
        """
        python {params.report_script} \
            {input} {output} \
            >{log} 2>&1
        """

# Create Hic File for Juicebox Visualisation
# Output file in MEDIUM format:
# <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
# Sort on both chromosome columns
rule clean_pairs:
    input:
        "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    output:
        temp("analysis_other/porec/{dataset}/pairs/{dataset}.pairs.for_juice")
    shell:
        """
        zcat {input} | \
        grep -v '^#' | \
        awk '{{OFS="\\t"; print $1,$6,$2,$3,0,$7,$4,$5,1,$11,$12}}' | \
        awk '{{OFS="\\t"; $2=($2=="+")?0:1; $6=($6=="+")?0:1; print}}' | \
        sort -k3,3V -k7,7V > {output}
        """

rule pairs_to_hic:
    input:
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.for_juice",
        chromsize = config['ref'] + ".chrom-size.txt"
    output:
        hic = "analysis_other/porec/{dataset}/hic/{dataset}.hic"
    params:
        resolutions = config.get("juicer_resolutions", "5000,10000,25000,50000,100000,500000")
    log:
        "logs/porec/juicer_tools_pre.{dataset}.log"
    shell:
        '''
        java -Xmx50G -jar /mnt/storage2/users/ahleucs1/tools/juicertools/juicer_tools_1.19.02.jar \
        pre {input.pairs} {output.hic} {input.chromsize} \
        -r {params.resolutions} > {log} 2>&1
        '''

# Create coolers for downstream analysis
rule pairs_to_cooler:
    input:
        fai = f"{config['ref']}.fai",
        pairs = "analysis_other/porec/{dataset}/pairs/{dataset}.pairs.gz"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.cool"
    wildcard_constraints:
        resolution = r"\d+"
    log:
        "logs/porec/pairs_to_cooler.{dataset}_{resolution}.log"
    conda:
        "../env/cooler.yml"
    threads: 2
    shell:
        """
        cooler cload pairs \
            -c1 2 -p1 3 -c2 4 -p2 5 \
            {input.fai}:{wildcards.resolution} \
            {input.pairs} \
            {output.cool} \
            >{log} 2>&1
        """

rule merge_mcools:
    input:
        expand("analysis_other/porec/{{dataset}}/cooler/{{dataset}}_{resolution}.cool", resolution=config.get("min_bin_width", "1000"))
    output:
        mcool = "analysis_other/porec/{dataset}/cooler/{dataset}.mcool"
    log:
        "logs/porec/merge_mcools.{dataset}.log"
    params:
        resolutions = config.get("cooler_resolutions", "1000,5000,10000,50000,100000"),
        prefix = lambda wc: wc.dataset
    conda:
        "../env/cooler.yml"
    threads: 2
    shell:
        """
        cooler zoomify \
            -r {params.resolutions} \
            -o {output.mcool} \
            {input}\
            >{log} 2>&1
        """

# Cooltools analysis
# 1. Cooler balance to correct matrix (ICE)
rule cooler_balance:
    input:
        "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.cool"
    output:
        "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_balanced.cool"
    log:
         "logs/porec/cooler_balance.{dataset}_{resolution}.log"
    threads:
        16
    conda:
        "../env/cooler.yml"
    shell:
        """
        cp {input} {output}
        cooler balance \
        --cis-only \
        --ignore-diags 2 \
        --mad-max 5 \
        --min-nnz 10 \
        --tol 1e-05 \
        --max-iters 500 \
        --nproc {threads} \
        {output} \
        >{log} 2>&1
        """   

#2. Compartments
rule eigs_cis:
    input:
        gc = config['gc_track'],
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_100000_balanced.cool"
    output:
        vecs="analysis_other/porec/{dataset}/compartments/{dataset}.cis.vecs.tsv",
        lam="analysis_other/porec/{dataset}/compartments/{dataset}.cis.lam.txt",
        bw="analysis_other/porec/{dataset}/compartments/{dataset}.cis.bw"
    params:
        prefix = "analysis_other/porec/{dataset}/compartments/{dataset}"
    log:
        "logs/porec/cooltools_eigs_cis.{dataset}.log"
    threads:
        2
    conda:
        "../env/cooltools.yml"
    shell:
        """
        cooltools eigs-cis \
        --phasing-track {input.gc}::GC \
        --clr-weight-name "weight" \
        --bigwig \
        --out-prefix {params.prefix} \
        {input.cool} \
        >{log} 2>&1
        """

#3. chromatin insulation
rule insulation_score:
    input: 
        cool =  "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_balanced.cool"
    output:
        insu = "analysis_other/porec/{dataset}/insulation/{dataset}_{resolution}.insulation.tsv",
        bw = "analysis_other/porec/{dataset}/insulation/{dataset}_{resolution}.insulation.tsv_bw_pseudo"        
    params:
        window = lambda wildcards: " ".join([str(mult * int(wildcards.resolution)) for mult in config['insu_window_multipliers']])
    log:
        "logs/porec/cooltools_insulation.{dataset}_{resolution}.log"
    conda:
        "../env/cooltools.yml"
    threads:
        8
    shell:
        """
        cooltools insulation \
        --threshold Li \
        --clr-weight-name "weight" \
        --output {output.insu} \
        --bigwig \
        --nproc {threads} \
        {input.cool} {params.window} \
        >{log} 2>&1

        touch {output.bw}
        """


# Does not currently do anything
rule hic_normalize:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.cool"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool"
    log:
        "logs/porec/hic_normalize.{dataset}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicNormalize \
            --matrices {input.cool} \
            --normalize smallest \
            -o {output.cool} \
            >{log} 2>&1
        """

rule hic_diagnostic_plot:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool"
    output:
        plot = "analysis_other/porec/{dataset}/qc/{dataset}_{resolution}_diagnostic.png"
    log:
        "logs/porec/hic_diagnostic.{dataset}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell: 
        """
        hicCorrectMatrix diagnostic_plot \
            --matrix {input.cool} \
            -o {output.plot} \
            >{log} 2>&1
        """

# new
rule hic_correct_matrix:
    input:
        cool =  "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_norm.cool",
        diagnostic = "logs/porec/hic_diagnostic.{dataset}.{resolution}.log"
    output:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        lower_file=temp("analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.lower"),
        upper_file=temp("analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}.upper")
    log:
        "logs/porec/hic_correct.{dataset}.{resolution}.log"
    params:
        correction_method = config.get("correction_method", "ICE")
        #filter_threshold = config.get("filter_threshold", "-1.5 4") 
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        echo "Start log ....." > {log}
        grep "mad threshold" {input.diagnostic} 2>> {log} | \
        sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g' 2>> {log} > {output.lower_file}
        echo -3*$(cat {output.lower_file}) | bc 2>>{log} > {output.upper_file}

        hicCorrectMatrix correct \
        --filterThreshold $(cat {output.lower_file}) $(cat {output.upper_file}) \
        --matrix {input.cool} \
        --correctionMethod {params.correction_method} \
        --outFileName {output.cool} \
        >>{log} 2>&1 || {{
            
            echo "Standard correction failed, trying more stringent filtering (lower z-score = -1)..." >>{log}
            
            # Try with more stringent filtering
            hicCorrectMatrix correct \
                --matrix {input.cool} \
                --correctionMethod {params.correction_method} \
                --outFileName {output.cool} \
                --filterThreshold -1.0 3.0 \
                --iterNum 200 \
                >>{log} 2>&1 || {{
                
                echo "Stringent correction failed, trying KR normalization..." >>{log}
                
                # Try Knight-Ruiz normalization as fallback
                hicCorrectMatrix correct \
                    --matrix {input.cool} \
                    --correctionMethod KR \
                    --outFileName {output.cool} \
                    --filterThreshold -3.0 3.0 \
                    >>{log} 2>&1 || {{
                    
                    echo "All correction methods failed, copying normalized matrix..." >>{log}
                    # As last resort, just copy the normalized matrix
                    cp {input.cool} {output.cool}
                }}
            }}
        }}
        
        """

# Plots
# 1. contact decay
rule hic_plot_dist_vs_counts:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        plot = "analysis_other/porec/{dataset}/qc/plot_vs_counts_{resolution}.png"
    log:
        "logs/porec/plot_dist_vs_counts.{dataset}.{resolution}.log"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotDistVsCounts \
            --matrices {input.cool} \
            -o {output.plot} \
            >{log} 2>&1
        """

# 2. whole chr +  eigs
rule fanc_plot_whole_chr_eigs_new:
    input:
        hic= "analysis_other/porec/{dataset}/hic/{dataset}_nonadj.hic",
        eigs="analysis_other/porec/{dataset}/compartments/{dataset}_nonadj.cis.bw"
    output:
        "analysis_other/porec/{dataset}/plots/whole_chr/{dataset}_{chrom}_None-None.png"
    params:
        chroms= lambda wc: wc.chrom,
        resolution=250000,
        name= lambda wc: wc.dataset,
        vmax=200
    conda:
        "../env/fanc.yml"
    log:
        "logs/porec/plotting/fanc_triangle_plot.whole_chr.{dataset}_{chrom}.log"
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

# 3. interesting regions + insulation (can add tad later)
# insulation resolution and chosen windows are hardcoded
# change if need be; possible values: 25000: 75000, 125000, 250000, 625000; 40000: 120000, 200000, 400000, 1000000
rule fanc_plot_imprinted_genes:
    input:
        insulation = "analysis_other/porec/{dataset}/insulation/{dataset}_25000.insulation.tsv_bw_pseudo", 
        hic1= lambda wc: f"analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1.hic", 
        hic2= lambda wc: f"analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/hic/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2.hic", 
    output:
        "analysis_other/porec/{dataset}/plots/imprinted_genes/{dataset}_{region}.png"
    params:
        insulation1= lambda wc: f"analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp1_25000.insulation.tsv.125000.bw", 
        insulation2= lambda wc: f"analysis_other/porec/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2/insulation/{re.sub(r'\.hp[12]', '', wc.dataset, count=1)}.hp2_25000.insulation.tsv.125000.bw",
        regions= lambda wc: wc.region.replace("_",":"),
        resolution=10000,
        name= lambda wc: wc.dataset,
        vmax=40
    conda:
        "../env/fanc.yml"
    log:
        "logs/porec/plotting/fanc_triangle_plot.imprinted_genes.{dataset}_{region}.log"
    shell:
        """
        fancplot \
        -n {params.name} \
        -o {output} {params.regions} \
        -p mirror \
        -uvmin 0 -uvmax {params.vmax} \
        -lvmin 0 -lvmax {params.vmax} \
        {input.hic1}@{params.resolution}@KR {input.hic2}@{params.resolution}@KR \
        -p line -l hp1 hp2 -c black red --fix-chromosome {params.insulation1} {params.insulation2} \
        >{log} 2>&1
        """

# architecture: tads and loops
rule hic_find_tads:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        tads = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_domains.bed",
        boundaries = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_boundaries.bed",
        score = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_score.bedgraph"
    log:
        "logs/porec/hic_find_tads.{dataset}.{resolution}.log"
    params:
        min_depth = lambda wc: max(config.get("tad_min_depth", 20000), 3*int(wc.resolution)),
        max_depth = config.get("tad_max_depth", 200000),
        step = lambda wc: max(config.get("tad_step", 10000), int(wc.resolution)),
        fdr_threshold = config.get("tad_fdr_threshold", 0.05),
        delta = config.get("tad_delta", 0.01),
        correction_factor_threshold = config.get("tad_correction_threshold", 1.5)
    conda:
        "../env/hicexplorer.yml"
    threads:
        4
    shell:
        """
        hicFindTADs \
            --matrix {input.cool} \
            --outPrefix analysis_other/porec/{wildcards.dataset}/tad/{wildcards.dataset}_{wildcards.resolution} \
            --minDepth {params.min_depth} \
            --maxDepth {params.max_depth} \
            --step {params.step} \
            --thresholdComparisons {params.fdr_threshold} \
            --delta {params.delta} \
            --correctForMultipleTesting fdr \
            --numberOfProcessors {threads} \
            >{log} 2>&1
        """

rule hic_plot_tads:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        tads = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_domains.bed"
    output:
        plot = "analysis_other/porec/{dataset}/tad/{dataset}_{resolution}_{region}_tads.png"
    log:
        "logs/porec/hic_plot_tads.{dataset}.{resolution}.{region}.log"
    params:
        region = lambda wc: wc.region,  # e.g., "chr1:1000000-5000000"
        dpi = config.get("plot_dpi", 300)
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicPlotMatrix \
            --matrix {input.cool} \
            --region {params.region} \
            --outFileName {output.plot} \
            --dpi {params.dpi} \
            --clearMaskedBins \
            >{log} 2>&1
        """

rule hic_detect_loops:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool"
    output:
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe",
    log:
        "logs/porec/hic_detect_loops.{dataset}.{resolution}.log"
    params:
        window_size = config.get("loop_window_size", 10),
        peak_width = config.get("loop_peak_width", 6),
        max_loop_distance = config.get("loop_max_distance", 2000000),
        pvalue_threshold = config.get("loop_pvalue_threshold", 0.05),
        fdr_threshold = config.get("loop_fdr_threshold", 0.1)
    threads: 4
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicDetectLoops \
            --matrix {input.cool} \
            --outFileName {output.loops} \
            --windowSize {params.window_size} \
            --peakWidth {params.peak_width} \
            --maxLoopDistance {params.max_loop_distance} \
            --pValuePreselection {params.pvalue_threshold} \
            --pValue {params.fdr_threshold} \
            --threads {threads} \
            >{log} 2>&1
        """

# Not used yet
rule hic_validate_viewpoints:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe"
    output:
        viewpoints = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_viewpoints.txt",
        interactions = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_interactions.txt"
    log:
        "logs/porec/hic_validate_viewpoints.{dataset}.{resolution}.log"
    params:
        resolution_val = lambda wc: wc.resolution,
        interaction_file = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_interactions.txt"
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicValidateLocations \
            --matrix {input.cool} \
            --locations {input.loops} \
            --method loops \
            --outFileName {output.viewpoints} \
            --outFileNameInteractions {params.interaction_file} \
            --resolution {params.resolution_val} \
            >{log} 2>&1
        """

# Not used yet
rule hic_aggregate_contacts:
    input:
        cool = "analysis_other/porec/{dataset}/cooler/{dataset}_{resolution}_corrected.cool",
        loops = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_loops.bedpe"
    output:
        aggregate_plot = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_aggregate_loops.png",
        aggregate_data = "analysis_other/porec/{dataset}/loops/{dataset}_{resolution}_aggregate_loops.npz"
    log:
        "logs/porec/hic_aggregate_contacts.{dataset}.{resolution}.log"
    params:
        range_around_loop = config.get("aggregate_range", 100000),
        number_of_bins = config.get("aggregate_bins", 30)
    conda:
        "../env/hicexplorer.yml"
    shell:
        """
        hicAggregateContacts \
            --matrix {input.cool} \
            --BED {input.loops} \
            --outFileName {output.aggregate_plot} \
            --outFilePrefixMatrix analysis_other/porec/{wildcards.dataset}/loops/{wildcards.dataset}_{wildcards.resolution}_aggregate_loops \
            --range {params.range_around_loop}:{params.range_around_loop} \
            --numberOfBins {params.number_of_bins} \
            --plotType 2d \
            >{log} 2>&1
        """

# Todo:
# Differential TAD: HP1 vs HP2
# Look at target regions:
#   - Imprinted genes