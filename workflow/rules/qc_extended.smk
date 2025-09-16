rule gqc:
    input:
        assembly = lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"],
        ref=config["ref_hg002_q100"]
    output:
        directory("analysis_other/GQC/{asm}")
    conda:
        "../env/gqc.yml"
    log:
        "logs/gqc/{asm}_gqc.log"
    threads: 10
    params:
        config = config['gqc']['config'],
        resourcedir = config['gqc']['resourcedir'],
        venv = config['gqc']['venv'],
        fastk_binfolder = config['gqc']['fastk_binfolder'],
    shell:
        """
        source {params.venv}/bin/activate
        export PATH=$PATH:{params.fastk_binfolder}
        GQC \
            --reffasta {input.ref} \
            --queryfasta {input.assembly} \
            --config {params.config} \
            -p {output} \
            -t {threads} \
            --assembly {wildcards.asm} \
            --benchmark HG002v1.1 \
            > {log} 2>&1
        """


rule repeatmasker:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        out="analysis_other/repeatmasker/{asm}/assembly.fasta.out",
        gff="analysis_other/repeatmasker/{asm}/assembly.fasta.out.gff",
    conda:
        "../env/repeatmasker.yml"
    log:
        "logs/repeatmasker/{asm}_repeatmasker.log"
    threads: 32
    benchmark:
        "runtimes/repeatmasker/{asm}.repeatmasker.txt"
    params:
        dfam_lib=config.get('dfam_db')
    shell:
        """
        RepeatMasker \
            -lib {params.dfam_lib} \
            -engine rmblast \
            -dir $(dirname {output.out}) \
            -pa {threads} \
            -gff \
            -q \
            -no_is \
            {input.fa} \
            > {log} 2>&1
        """


rule analyze_repeatmasker:
    input:
        rm_out="analysis_other/repeatmasker/{asm}/assembly.fasta.out"
    output:
        summary="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_sequence_summary.csv",
        bed = "analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucflag.bed"
    conda:
        "../env/R.yml"
    log:
        "logs/analyze_repeatmasker/{asm}_analyze_repeatmasker.log"
    shell:
        """
        # Run the R script to analyze the RepeatMasker output
        Rscript workflow/scripts/16_analyze_repeatmasker.R \
            --input {input.rm_out} \
            --output $(dirname {output.summary})/{wildcards.asm}   \
            > {log} 2>&1
        """

rule nucflag:
    input:
        bam="data/mapped/{asm}/{asm}.HQ.asm.bam",
        rm_bed="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucflag.bed"
    output:
        status="analysis_other/nucflag/{asm}/nucflag_status.bed",
        misasm="analysis_other/nucflag/{asm}/nucflag_misasm.bed",
        plots=directory("analysis_other/nucflag/{asm}/plots")
    conda:
        "../env/nucflag.yml"
    log:
        "logs/nucflag/{asm}_nucflag.log"
    threads: 12
    shell:
        """
        nucflag \
            -i {input.bam} \
            -d {output.plots} \
            --overlay_regions {input.rm_bed} \
            --threads {threads} \
            --output_status {output.status} \
            --output_misasm {output.misasm} \
            > {log} 2>&1

        # Replace colons with underscores in plot filenames
        for file in {output.plots}/*; do
            if [[ -f "$file" ]]; then
                newname=$(echo "$file" | tr ':' '_')
                if [[ "$file" != "$newname" ]]; then
                    mv "$file" "$newname"
                fi
            fi
        done
        """

rule asm_index:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        "analysis_other/flagger/{asm}/assembly.index.bed"
    log:
        "logs/flagger/{asm}/asm_index.log"
    threads:
        1
    shell:
        """
        samtools faidx {input.fa} && \
        awk 'BEGIN{{OFS="\t"}} {{print $1, 0, $2}}' {input.fa}.fai > {output} 2>{log}
        """

rule create_index_json:
    input:
        "analysis_other/flagger/{asm}/assembly.index.bed"
    output:
        "analysis_other/flagger/{asm}/assembly.index.json"
    log:
        "logs/flagger/{asm}/create_json.log"
    params:
        workdir = workflow.basedir
    threads:
        1
    shell:
        """
        echo "{{" > {output}
        echo \\"whole_genome\\" : \\"/mnt/{input}\\" >> {output}
        echo "}}" >> {output}
        """

rule flagger_create_cov:
    input:
        bam = "data/mapped/{asm}/{asm}.HQ.asm.bam",
        idx = "analysis_other/flagger/{asm}/assembly.index.json"
    output:
        cov = "analysis_other/flagger/{asm}/coverage_file.cov.gz",
    log:
        "logs/flagger/{asm}/create_cov.log"
    threads:
        16
    params:
        bam2cov = config['bam2cov'], 
        workdir = workflow.basedir
    shell:
        """
        docker run --rm \
            --user $(id -u):$(id -g) \
            -v$(dirname {params.workdir}):/mnt \
            {params.bam2cov} \
                --bam /mnt/{input.bam} \
                --output /mnt/{output.cov} \
                --annotationJson /mnt/{input.idx} \
                --threads {threads} \
                --baselineAnnotation "whole_genome" \
            >{log} 2>&1
        """

rule flagger:
    input:
        cov = "analysis_other/flagger/{asm}/coverage_file.cov.gz"
    output:
        tsv = "analysis_other/flagger/{asm}/prediction_summary_final.tsv",
        bed = "analysis_other/flagger/{asm}/final_flagger_prediction.bed"
    log:
        "logs/flagger/{asm}/hmm_flager.log"
    threads:
        16
    params:
        hmm_flagger = config['hmm_flagger'],
        alpha_tsv = config['hmm_flagger_alpha_tsv'],
        workdir = workflow.basedir
    shell:
        """
        docker run --rm \
            --user $(id -u):$(id -g) \
            -v$(dirname {params.workdir}):/mnt \
            {params.hmm_flagger} \
                --input /mnt/{input.cov} \
                --outputDir /mnt/$(dirname {output.tsv}) \
                --trackName {wildcards.asm}_flagger_prediction \
                --threads {threads} \
        >{log} 2>&1
        """

rule create_plot:
    input:
        contig_file="analysis_other/assembly_info/{asm}/contigs.txt",
        fai_file="analysis_other/assembly_info/{asm}/assembly.fasta.fai",
        feature_bed="analysis_other/features/{asm}/features.bed"
    output:
        plot="analysis_other/plots/{asm}/{asm}_genomic_features.pdf"
    conda:
        "../env/R.yml"
    log:
        "logs/create_plot/{asm}_create_plot.log"
    shell:
        """
        mkdir -p $(dirname {output.plot})
        
        Rscript workflow/scripts/20_plot_genomic_features.R \
            --contig_file {input.contig_file} \
            --fai_file {input.fai_file} \
            --feature_files {input.feature_bed} \
            --output_prefix analysis_other/plots/{wildcards.asm}/{wildcards.asm} \
            > {log} 2>&1
        """