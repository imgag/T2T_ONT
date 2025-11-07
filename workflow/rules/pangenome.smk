rule all_pangenome:
    input:
        "analysis_other/pangenome/multiqc_report.html"

rule filter_assembly:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        "analysis_other/pangenome/filtered/{asm}.fasta"
    log:
        "logs/filter_assembly/{asm}.log"
    conda:
        "../env/samtools.yml"
    shell:
        """
        # Create temporary files for processing
        temp_fasta=$(mktemp)
        temp_lengths=$(mktemp)
        
        # First, calculate sequence lengths
        echo "Calculating contig lengths..." > {log}
        awk '
        /^>/ {{
            if (seq_name != "") {{
                print seq_name "\\t" length(sequence)
            }}
            seq_name = substr($0, 2)
            sequence = ""
        }}
        !/^>/ {{
            sequence = sequence $0
        }}
        END {{
            if (seq_name != "") {{
                print seq_name "\\t" length(sequence)
            }}
        }}
        ' {input.fa} > $temp_lengths 2>>{log}
        
        # Log original statistics
        total_contigs=$(wc -l < $temp_lengths)
        echo "Original contigs: $total_contigs" >> {log}
        
        # Count unassigned contigs
        unassigned_count=$(grep "^unassigned" $temp_lengths | wc -l)
        echo "Unassigned contigs: $unassigned_count" >> {log}
        
        # Count short contigs
        short_count=$(awk '$2 < 250000' $temp_lengths | wc -l)
        echo "Short contigs (<250kb): $short_count" >> {log}
        
        # Count contigs to keep
        keep_count=$(awk '$2 >= 250000 && !match($1, /^unassigned/)' $temp_lengths | wc -l) 
        echo "Contigs to keep: $keep_count" >> {log}
        
        # Filter the FASTA file
        echo "Filtering FASTA file..." >> {log}
        awk -v lengths_file="$temp_lengths" '
        BEGIN {{
            # Load contig lengths and filter criteria
            while ((getline line < lengths_file) > 0) {{
                split(line, parts, "\\t")
                contig_name = parts[1]
                contig_len = parts[2]
                
                # Keep contigs that are >= 250kb and not unassigned
                if (contig_len >= 250000 && !match(contig_name, /^unassigned/)) {{
                    keep[contig_name] = 1
                }}
            }}
            close(lengths_file)
        }}
        /^>/ {{
            seq_name = substr($0, 2)
            if (seq_name in keep) {{
                print_seq = 1
                print $0
            }} else {{
                print_seq = 0
            }}
            next
        }}
        print_seq == 1 {{ print }}
        ' {input.fa} > {output} 2>>{log}
        
        # Log final statistics
        final_contigs=$(grep "^>" {output} | wc -l)
        echo "Final contigs: $final_contigs" >> {log}
        
        # Calculate total lengths
        original_bp=$(awk '{{sum += $2}} END {{print sum}}' $temp_lengths)
        final_bp=$(awk '/^>/ {{if (seq_name != "") {{total += length(sequence)}} seq_name = substr($0, 2); sequence = ""}} !/^>/ {{sequence = sequence $0}} END {{if (seq_name != "") {{total += length(sequence)}} print total}}' {output})
        
        echo "Original total length: $original_bp bp" >> {log}
        echo "Final total length: $final_bp bp" >> {log}
        
        # Cleanup
        rm -f $temp_fasta $temp_lengths
        """

rule rename_assembly_header:
# PanSN naming convention (https://github.com/pangenome/PanSN-spec)
# [sample_name][delim][haplotype_id][delim][contig_id]
    input:
        fa="analysis_other/pangenome/filtered/{asm}.fasta"
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