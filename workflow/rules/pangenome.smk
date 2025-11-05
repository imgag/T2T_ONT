rule all_pangenome:
    input:
        "analysis_other/pangenome/multiqc_report.html"

rule rename_assembly_header:
# PanSN naming convention (https://github.com/pangenome/PanSN-spec)
# [sample_name][delim][haplotype_id][delim][contig_id]
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        "analysis_other/pangenome/assemblies/{asm}.fasta",
    log:
        "logs/rename_assembly/{asm}.log"
    shell:
        """
        awk -v sample="{wildcards.asm}" '
        /^>/ {{
            # Extract original header without ">"
            header = substr($0, 2)
            
            # Split on hyphen to get haplotype and contig number
            # Example: haplotype1-0000001 -> haplotype1, 0000001
            split(header, parts, "-")
            
            # Extract haplotype number (1 or 2)
            hap = substr(parts[1], length(parts[1]), 1)
            
            # Get contig number
            contig = parts[2]
            
            # Create PanSN format: sample#hap#contig
            print ">" sample "#" hap "#" contig
            next
        }}
        {{ print }}
        ' {input.fa} > {output} 2>{log}
        """

rule merge_assemblies:
    input:
        expand("analysis_other/pangenome/assemblies/{asm}.fasta", asm=finished_samples),
    output:
        "analysis_other/pangenome/all_assemblies.fasta.bgz"
    log:
        "logs/pangenome/merge_assemblies.log"
    threads: 8
    conda:
        "../env/samtools.yml"
    shell:
        """
        cat {input} | bgzip -@{threads} -c > {output} 2>{log}
        samtools faidx {output}
        """

rule nfcore_pangenome:
    input:
        assemblies = "analysis_other/pangenome/all_assemblies.fasta.bgz",
    output:
        report = "analysis_other/pangenome/multiqc_report.html",
    threads: 64
    benchmark:
        "runtimes/pangenome/nfcore_pangenome.txt"
    log:
        "logs/pangenome/nfcore_pangenome.log"
    params:
        nf_plugins_dir = "$PWD/.nextflow/plugins"
    shell:
        """
        tmp_config=$(mktemp)
        echo "process {{ resourceLimits = [cpus : {threads}] }}" > $tmp_config

        export NXF_PLUGINS_DIR={params.nf_plugins_dir}
        .bin/nextflow run bin/nf-core-pangenome_1.1.3/1_1_3 \
            --input {input.assemblies} \
            --n_haplotypes 2 \
            --outdir analysis_other/pangenome
            -profile singularity \
            -c $tmp_config \
            > {log} 2>&1
        """