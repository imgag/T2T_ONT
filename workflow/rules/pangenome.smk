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

rule rename_hg002_header:
# PanSN naming convention for HG002 reference
# Convert chr1_MATERNAL -> HG002#1#chr1, chr1_PATERNAL -> HG002#2#chr1
    input:
        fa=config['ref_hg002_q100']  # Adjust path as needed
    output:
        "analysis_other/pangenome/assemblies/HG002.fasta"
    log:
        "logs/rename_assembly/HG002.log"
    shell:
        """
        zcat {input.fa} | awk -v sample="HG002" '
        /^>/ {{
            # Extract original header without ">"
            header = substr($0, 2)
            
            # Split on underscore to get chromosome and haplotype
            # Example: chr1_MATERNAL -> chr1, MATERNAL
            split(header, parts, "_")
            
            chr = parts[1]
            haplotype = parts[2]
            
            # Map MATERNAL to haplotype 1, PATERNAL to haplotype 2
            if (haplotype == "MATERNAL") {{
                hap = "1"
            }} else if (haplotype == "PATERNAL") {{
                hap = "2"
            }} else {{
                # If no haplotype specified (e.g., unphased contigs), default to 0
                hap = "0"
            }}
            
            # Create PanSN format: sample#hap#chr
            print ">" sample "#" hap "#" chr
            next
        }}
        {{ print }}
        ' > {output} 2>{log}
        """

rule merge_assemblies:
    input:
        expand("analysis_other/pangenome/assemblies/{asm}.fasta", asm=finished_samples) +
        ["analysis_other/pangenome/assemblies/HG002.fasta"],
    output:
        "analysis_other/pangenome/all_assemblies.fasta.gz"
    log:
        "logs/pangenome/merge_assemblies.log"
    threads: 8
    conda:
        "../env/samtools.yml"
    shell:
        """
        cat {input} | bgzip -@{threads} -c > {output} 2>{log}
        samtools faidx {output} 2>>{log}
        """

rule nfcore_pangenome:
    input:
        assemblies = "analysis_other/pangenome/all_assemblies.fasta.gz",
    output:
        report = "analysis_other/pangenome/multiqc_report.html",
    threads: 64
    benchmark:
        "runtimes/pangenome/nfcore_pangenome.txt"
    log:
        "logs/pangenome/nfcore_pangenome.log"
    params:
        nf_plugins_dir = "$PWD/.nextflow/plugins",
        config  = "bin/nf-core-pangenome_1.1.3/1_1_3/imgag.nextflow.config",
    shell:
        """

        export NXF_PLUGINS_DIR={params.nf_plugins_dir}
        export NXF_OFFLINE='true'
        export NXF_SINGULARITY_CACHEDIR="bin/nf-core-pangenome_1.1.3/singularity-images"

        ./bin/nextflow run bin/nf-core-pangenome_1.1.3/1_1_3 \
            --input {input.assemblies} \
            --n_haplotypes 2 \
            --outdir analysis_other/pangenome \
            -profile singularity \
            -c {params.config} \
            > {log} 2>&1
        """