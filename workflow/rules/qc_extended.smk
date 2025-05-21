rule gqc:
    input:
        ref=config["ref"],
        asm=get_assembly_output,
    output:
        directory("analysis_other/GQC/{asm}")
    conda:
        "../env/gqc.yml"
    log:
        "logs/gqc/{isphased}_{tool}_{asm}.log"
    threads: 10
    params:
        config = config['gqc']['config'],
        resourcedir = config['gqc']['resourcedir'],
        venv = config['gqc']['venv'],
        fastk_binfolder = config['gqc']['fastk_binfolder'],
        ref = config['ref'],
    shell:
        """
        source {params.venv}/bin/activate
        export PATH=$PATH:{params.fastk_binfolder}
        GQC \
            --reffasta {input.ref} \
            --queryfasta {input.asm} \
            --config {params.config} \
            -p {output} \
            -t {threads} \
            --assembly {wildcards.asm} \
            --benchmark HG002v1.1 \
            > {log} 2>&1
        """

rule map_ul_to_asm:
    input:
        fa=get_assembly_output
    output
        bam=
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_ul_to_asm/{asm}.log"
    threads: 40
    shell:
        """
        
        """


rule repeatmasker:
    input:
        fa=get_assembly_output
    output:
        out="analysis_other/repeatmasker/{asm}/assembly.fasta.out",
        gff="analysis_other/repeatmasker/{asm}/assembly.fasta.out.gff",
        html="analysis_other/repeatmasker/{asm}/assembly.fasta.html"
    conda:
        "../env/repeatmasker.yml"
    log:
        "logs/repeatmasker/{isphased}_{tool}_{asm}.log"
    threads: 64
    shell:
        """
        # Run RepeatMasker
        RepeatMasker \
            -species human \
            -dir $(dirname {output.out}) \
            -pa {threads} \
            -gff \
            -html \
            {input.fa} \
            > {log} 2>&1
        """

rule analyze_repeatmasker:
    input:
        rm_out="analysis_other/repeatmasker/{asm}/assembly.fasta.out"
    output:
        summary="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_sequence_symmary.csv",
        bed = "analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucplod.bed"
    conda:
        "../env/R.yml"
    log:
        "logs/analyze_repeatmasker/{isphased}_{tool}_{asm}.log"
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
        bam="data/bam/{asm}.HQ_herro.50x.bam",
        rm_bed="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucplod.bed"
    output:
        status="analysis_other/nucflag/{asm}/nucflag_status.bed",
        misasm="analysis_other/nucflag/{asm}/nucflag_misasm.bed",
        plots=directory("analysis_other/nucflag/{asm}/plots")
    conda:
        "../env/nucflag.yml"
    log:
        "logs/nucflag/{isphased}_{tool}_{asm}.log"
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