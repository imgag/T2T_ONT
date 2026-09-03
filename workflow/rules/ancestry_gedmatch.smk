rule all_ancestry_gedmatch:
    input:
        expand("analysis_other/gedmatch/{asm}/{asm}_gedmatch_upload.txt", asm=finished_samples)

rule extract_array_snps_from_manifest:
    #Extract SNP list from Illumina Global Screening Array manifest CSV
    #Include rsID, chromosome, position, and REF alleles
    input:
        manifest_csv=config["array_manifest_csv"]
    output:
        rsids="data/dbsnp/array_snps.txt",
        positions="data/dbsnp/array_snp_positions_with_ref.txt"
    log:    
        "logs/ancestry_gedmatch/extract_array_snps_from_manifest.log"
    shell:
        """
        echo "Extracting SNPs from Illumina array manifest" > {log}

        # Skip first 8 lines (header + column names), extract Name column (field 2)
        # Filter for rsIDs only
        tail -n +9 {input.manifest_csv} | cut -d',' -f2 | grep '^rs' | sort -u > {output.rsids}

        # Extract rsID, chromosome (field 10), position (field 11), and SNP (field 4)
        # SNP field format is [REF/ALT] - extract REF allele
        tail -n +9 {input.manifest_csv} | \
            awk -F',' '$2 ~ /^rs/ {{
                rsid = $2
                chr = $10
                pos = $11
                snp = $4
                # Extract REF allele from [A/G] format - take first allele
                gsub(/[\[\]]/, "", snp)
                split(snp, alleles, "/")
                ref = alleles[1]
                print rsid "\\t" chr "\\t" pos "\\t" ref
            }}' | \
            grep -v '^rs.*\\t\\t' > {output.positions} 2>> {log}

        echo "Array SNPs extracted: $(wc -l < {output.rsids})" >> {log}
        echo "Array positions with REF: $(wc -l < {output.positions})" >> {log}
        """



rule dipcall_hg19:
    # Small Variants and Indels with dipcall against hg19 ref
    input:
        pat_fa="assembly/output/verkko/{asm}/assembly.haplotype1.fasta",
        mat_fa="assembly/output/verkko/{asm}/assembly.haplotype2.fasta",
        ref_fa=config['ref_grch37'],
        sex="assembly/qc/phased_verkko/{asm}/sample_sex.txt",
    output:
        makefile="assembly/variants/hg19/{asm}/small_variants.dip.mak",
        vcf="assembly/variants/hg19/{asm}/small_variants.dip.vcf.gz",
    params:
        sex=lambda wc: "-x " + config["X_PAR_file"] if open(f"assembly/qc/phased_verkko/{wc.asm}/sample_sex.txt").read().strip() == "male" else "",
        run_dipcall=config["run-dipcall"],
    log:
        "logs/dipcall/{asm}.log",
    threads: 8
    shell:
        """
        {params.run_dipcall} \
            $(dirname {output.vcf})/small_variants \
            {input.ref_fa} \
            {input.pat_fa} {input.mat_fa} \
            {params.sex} -t {threads} \
            > {output.makefile} 2>{log}
        make -j2 -f {output.makefile} >{log} 2>&1
        tabix {output.vcf}
        """

rule fix_dbsnp_chromosome_names:
    """
    Convert RefSeq chromosome names (NC_000001.10) to standard format (chr1)
    using NCBI assembly report
    """
    input:
        vcf=config['dbsnp_hg19'],
        assembly_report=config['dbsnp_hg19_report']
    output:
        vcf="data/dbsnp/hg19_dbsnp_fixed.vcf.gz",
        tbi="data/dbsnp/hg19_dbsnp_fixed.vcf.gz.tbi"
    log:
        "logs/ancestry_gedmatch/fix_dbsnp_chromosome_names.log"
    conda:
        "../env/bcftools.yml"
    shell:
        """
        echo "Fixing chromosome names in dbSNP VCF" > {log}
        
        # Extract RefSeq to UCSC mapping from assembly report
        # Column 7 = RefSeq-Accn, Column 10 = UCSC-style-name
        grep -v "^#" {input.assembly_report} | \
            awk -F'\\t' '$7 != "na" && $7 != "" && $10 != "na" && $10 != "" {{print $7"\\t"$10}}' \
            > /tmp/chr_rename.txt 2>> {log}
        
        echo "Chromosome mapping created:" >> {log}
        head -10 /tmp/chr_rename.txt >> {log}
        echo "Total mappings: $(wc -l < /tmp/chr_rename.txt)" >> {log}
        
        chrom_list=$(cut -f1 /tmp/chr_rename.txt | tr '\\n' ',' | sed 's/,$//')
        
        echo "Filtering and renaming chromosomes..." >> {log}
        
        # Filter to only chromosomes with mappings, then rename and sort
        bcftools view \
            --targets "$chrom_list" \
            {input.vcf} 2>> {log} \
        | grep -v "^##contig=<ID=na>" 2>> {log} \
        | bcftools annotate \
            --rename-chrs /tmp/chr_rename.txt \
            --output-type v \
            - 2>> {log} \
        | bcftools sort \
            --output-type z \
            --output {output.vcf} \
            - 2>> {log}
        
        tabix -p vcf {output.vcf} 2>> {log}
        """

rule annotate_rsids_hg19:
    # Annotate dipcall hg19 VCF with rsIDs from dbSNP VCF 
    input:
        vcf="assembly/variants/hg19/{asm}/small_variants.dip.vcf.gz",
        dbsnp_vcf="data/dbsnp/hg19_dbsnp_fixed.vcf.gz",
    output:
        vcf="analysis_other/gedmatch/{asm}/variants_hg19_annotated.vcf.gz",
        tbi="analysis_other/gedmatch/{asm}/variants_hg19_annotated.vcf.gz.tbi"
    log:
        "logs/ancestry_gedmatch/{asm}/annotate_rsids_hg19.log"
    conda:
        "../env/bcftools.yml"
    threads: 4
    shell:
        """
        echo "Starting rsID annotation for {wildcards.asm}" > {log}
        
        bcftools annotate \
            --annotations {input.dbsnp_vcf} \
            --columns ID \
            --output-type z \
            --threads {threads} \
            {input.vcf} > {output.vcf} 2>> {log}
        
        tabix -p vcf {output.vcf} 2>> {log}
        
        # Log statistics
        echo "Total variants: $(bcftools view -H {output.vcf} | wc -l)" >> {log}
        echo "Variants with rsID: $(bcftools view -H {output.vcf} | awk '$3 ~ /^rs/' | wc -l)" >> {log}
        """


rule create_gedmatch_template:
    """
    Create template Gedmatch upload file with all array positions as REF calls
    """
    input:
        positions="data/dbsnp/array_snp_positions_with_ref.txt"
    output:
        "data/dbsnp/gedmatch_template.txt"
    log:
        "logs/ancestry_gedmatch/create_gedmatch_template.log"
    shell:
        """
        echo "Creating Gedmatch template with REF calls" > {log}
        
        # Write header
        echo -e "rsid\\tchromosome\\tposition\\tallele1\\tallele2" > {output}
        
        # Add all positions with REF/REF genotype
        awk '{{print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$4}}' {input.positions} >> {output}
        
        total=$(tail -n +2 {output} | wc -l)
        echo "Template created with $total REF/REF calls" >> {log}
        """

rule convert_to_gedmatch:
    #Update template with actual genotypes from VCF
    #Replaces REF/REF calls with actual variants where present
    input:
        vcf="analysis_other/gedmatch/{asm}/variants_hg19_annotated.vcf.gz",
        template="data/dbsnp/gedmatch_template.txt"
    output:
        "analysis_other/gedmatch/{asm}/{asm}_gedmatch_upload.txt"
    log:
        "logs/ancestry_gedmatch/{asm}/convert_to_gedmatch.log"
    conda:
        "../env/hapdiff.yml"
    script:
        "../scripts/43_vcf_to_gedmatch.py"