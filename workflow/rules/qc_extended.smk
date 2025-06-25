rule gqc:
    input:
        ref=config["ref_hg002_q100"],
        query=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"}),
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
            --queryfasta {input.query} \
            --config {params.config} \
            -p {output} \
            -t {threads} \
            --assembly {wildcards.asm} \
            --benchmark HG002v1.1 \
            > {log} 2>&1
        """

#rule map_ul_to_asm:
#    input:
#        fa=get_assembly_output
#    output
#        bam=
#    conda:
#        "../env/minimap2.yml"
#    log:
#        "logs/map_ul_to_asm/{asm}.log"
#    threads: 40
#    shell:
#        """
#        
#        """

rule repeatmasker:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})
    output:
        out="analysis_other/repeatmasker/{asm}/assembly.fasta.out",
        gff="analysis_other/repeatmasker/{asm}/assembly.fasta.out.gff",
        html="analysis_other/repeatmasker/{asm}/assembly.fasta.html"
    conda:
        "../env/repeatmasker.yml"
    log:
        "logs/repeatmasker/{asm}_repeatmasker.log"
    threads: 64
    benchmark:
        "runtimes/repeatmasker/{asm}.repeatmasker.txt"
    shell:
        """
        # Run RepeatMasker
        RepeatMasker \
            -species human \
            -dir $(dirname {output.out}) \
            -pa {threads} \
            -gff \
            {input.fa} \
            > {log} 2>&1
        """

rule analyze_repeatmasker:
    input:
        rm_out="analysis_other/repeatmasker/{asm}/assembly.fasta.out"
    output:
        summary="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_sequence_summary.csv",
        bed = "analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucplot.bed"
    conda:
        "../env/R.yml"
    log:
        "logs/analyze_repeatmasker/{asm}_analyze_repeatmasker.log"
    shell:
        """
        # Run the R script to analyze the RepeatMasker output
        Rscript workflow/scripts/16_analyze_repeatmasker.R \
            --input {input.rm_out} \
            --output {output.summary} \
            > {log} 2>&1
        """

rule nucflag:
    input:
        bam=lambda wc: get_assembly_input(wc).get("ul"),
        rm_bed="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucplot.bed"
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
