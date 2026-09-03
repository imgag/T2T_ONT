#ANCESTRY_TOOLS = ["iadmix"]

#LOCAL_ANCESTRY_TOOLS = ["gnomix", "rfmix"]
LOCAL_ANCESTRY_TOOLS = ["rfmix"]

#THOUSAND_G_POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

rule all_ancestry:
    input:
        # Merge sample and reference VCFs
        "analysis_other/ancestry/processed_vcf/merged_samples_ref.vcf.gz",
        # Calculate PCA with plink
        "analysis_other/ancestry/pca/merge_ref_samples.eigenvec",
        expand("analysis_other/ancestry/pca/pca_plot.{population}.pdf", population=["POP", "SUPERPOP"]),
        # Global ancestry
        #"analysis_other/ancestry/global/admixture/results.done",        
        # Local ancestry results  
        expand("analysis_other/ancestry/local/{tool}/results.done", tool=LOCAL_ANCESTRY_TOOLS),
        # RFMIX tagore plots
        "analysis_other/ancestry/local/rfmix/tagore_plots.done"


# ------------------ #
# Data Preprocessing #
# ------------------ #


# Filter and normalize single sample VCF
rule process_single_sample_vcf:
    input:
        vcf = "assembly/variants/{sample}/phased_verkko/small_variants.vcf.gz",
        ref = config["ref"]
    output:
        "analysis_other/ancestry/processed_vcf/samples/{sample}_phased.vcf.gz"
    conda:
        "../env/bcftools.yml" 
    log:
        "logs/ancestry/filter_vcf/process_single_sample_vcf_{sample}.log"
    shell:
        """
        # Normalize and index the single sample VCF
        bcftools norm \
            --fasta-ref  {input.ref} \
            --multiallelics -both \
            --check-ref w \
            --output-type u \
            {input.vcf} 2>{log} | \
        bcftools view \
            --types snps,indels \
            --exclude 'ALT~"\\*" || (strlen(REF)>50 || strlen(ALT)>50)' \
            --min-alleles 2 \
            --max-alleles 2 \
            --output-type z \
            --output {output} \
            >>{log} 2>&1
        
        tabix -p vcf {output} 2>>{log}
        """


# Extract ref sites (CHR, POS, REF, ALT) from 1000G VCF
rule extract_ref_sites:
    input:
        ref=config['1000G_phased_ref_vcf_biallelic']
    output:
        ref_sites="analysis_other/ancestry/processed_vcf/ref_sites.tsv"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/extract_ref_sites.log"
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.ref} > {output.ref_sites}
        """

# Subset each sample to reference sites only (one-sided join)
rule subset_samples_to_ref_sites:
    input:
        vcf="analysis_other/ancestry/processed_vcf/samples/{sample}_phased.vcf.gz",
        ref_sites="analysis_other/ancestry/processed_vcf/ref_sites.tsv"
    output:
        vcf="analysis_other/ancestry/processed_vcf/tmp/{sample}_phased.ref_sites.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/subset_to_ref_sites_{sample}.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/processed_vcf/tmp
        
        bcftools view \
            -T {input.ref_sites} \
            -Oz -o {output.vcf} \
            {input.vcf} >{log} 2>&1
        
        tabix -p vcf {output.vcf} >>{log} 2>&1
        """

# Build merge list (subsetted samples + reference)
rule build_merge_list:
    input:
        vcfs=lambda wc: expand(
            "analysis_other/ancestry/processed_vcf/tmp/{sample}_phased.ref_sites.vcf.gz",
            sample=[s for s in finished_samples if s in asm_samples]
        ),
        ref=config['1000G_phased_ref_vcf_biallelic']
    output:
        list="analysis_other/ancestry/processed_vcf/merge_list.txt"
    log:
        "logs/ancestry/build_merge_list.log"
    run:
        with open(output.list, "w") as f:
            for v in input.vcfs:
                f.write(v + "\n")
            f.write(input.ref + "\n")

# Create merged VCF of all samples + reference for PCA and iAdmix
rule merge_ref_and_samples:
    input:
        list="analysis_other/ancestry/processed_vcf/merge_list.txt"
    output:
        merged="analysis_other/ancestry/processed_vcf/merged_samples_ref.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/sample_merge_and_ref_vcfs.log"
    shell:
        """
        bcftools merge \
            --merge none \
            --missing-to-ref \
            --file-list {input.list} \
            --output-type z \
            --output {output.merged} \
            >{log} 2>&1
        
        tabix -p vcf {output.merged} >>{log} 2>&1
        """

# ----- #
#  PCA  #
# ----- # 

# Create PSAM file from metadata table and full merged VCF
rule create_psam_file:
    input:
        vcf="analysis_other/ancestry/processed_vcf/merged_samples_ref.vcf.gz",
        panel=config['1000G_metadata_file']
    output:
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/create_psam_file.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/plink/reference
        
        # Create PSAM file from merged VCF and metadata
        python3 workflow/scripts/26_convert_metadata_to_psam.py \
            --vcf {input.vcf} \
            --panel {input.panel} \
            --output {output.psam} \
            >{log} 2>&1
        """

# Convert VCF to Plink 2 format (PGEN) for PCA
rule vcf_to_plink:
    input:
        vcf="analysis_other/ancestry/processed_vcf/merged_samples_ref.vcf.gz",
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam"
    output:
        pgen="analysis_other/ancestry/plink/merge_ref_samples.pgen",
        pvar="analysis_other/ancestry/plink/merge_ref_samples.pvar",
        psam="analysis_other/ancestry/plink/merge_ref_samples.psam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/samples_to_plink.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/plink
        
        plink2 \
            --vcf {input.vcf} \
            --psam {input.psam} \
            --make-pgen \
            --out analysis_other/ancestry/plink/merge_ref_samples \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --set-missing-var-ids @:# \
            --threads {threads} \
            >{log} 2>&1
        """


# Filter datasets with QC filters
rule qc_filter_plink:
    input:
        pgen="analysis_other/ancestry/plink/{prefix}.pgen",
        pvar="analysis_other/ancestry/plink/{prefix}.pvar",
        psam="analysis_other/ancestry/plink/{prefix}.psam"
    output:
        pgen="analysis_other/ancestry/plink/{prefix}.filtered.pgen",
        pvar="analysis_other/ancestry/plink/{prefix}.filtered.pvar",
        psam="analysis_other/ancestry/plink/{prefix}.filtered.psam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/qc_filter_{prefix}.log"
    threads: 32

    # Unused filter:         
    # MIND_THRESHOLD=0.6     # Remove samples with >60% missing data (very lenient)
    # --mind $MIND_THRESHOLD 

    shell:
        """
        # QC filters - lenient for ancestry analysis
        MAF_THRESHOLD=0.01     # Keep rare variants for ancestry
        GENO_THRESHOLD=0.15    # Lenient missing rate
        HWE_THRESHOLD=1e-3     # Lenient HWE for population structure with families

        echo "Applying QC filters to {wildcards.prefix}..." >{log}
        
        plink2 \
            --pfile analysis_other/ancestry/plink/{wildcards.prefix} \
            --maf $MAF_THRESHOLD \
            --geno $GENO_THRESHOLD \
            --hwe $HWE_THRESHOLD \
            --make-pgen \
            --out analysis_other/ancestry/plink/{wildcards.prefix}.filtered \
            --threads {threads} \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 1000 \
            --sort-vars \
            >{log} 2>&1
        
        # Report filtering statistics
        echo "Filtering completed for {wildcards.prefix}" >>{log}
        echo "Variants after filtering: $(wc -l < analysis_other/ancestry/plink/{wildcards.prefix}.filtered.pvar | tail -n +2)" >>{log}
        echo "Samples after filtering: $(wc -l < analysis_other/ancestry/plink/{wildcards.prefix}.filtered.psam | tail -n +2)" >>{log}
        """

# PCA for population structure using PGEN format
rule run_pca:
    input:
        pgen="analysis_other/ancestry/plink/merge_ref_samples.filtered.pgen",
        pvar="analysis_other/ancestry/plink/merge_ref_samples.filtered.pvar", 
        psam="analysis_other/ancestry/plink/merge_ref_samples.filtered.psam"
    output:
        eigenval="analysis_other/ancestry/pca/merge_ref_samples.eigenval",
        eigenvec="analysis_other/ancestry/pca/merge_ref_samples.eigenvec"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/run_pca.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/pca
        
        # LD pruning for PCA
        plink2 \
            --pfile analysis_other/ancestry/plink/merge_ref_samples.filtered \
            --indep-pairwise 50 10 0.2 \
            --out analysis_other/ancestry/pca/ld_pruned \
            --threads {threads} \
            >{log} 2>&1
        
        # Run PCA
        plink2 \
            --pfile analysis_other/ancestry/plink/merge_ref_samples.filtered \
            --extract analysis_other/ancestry/pca/ld_pruned.prune.in \
            --pca 20 \
            --out analysis_other/ancestry/pca/merge_ref_samples \
            --threads {threads} \
            >{log} 2>&1
        """


rule plot_pca:
    input:
        matrix="analysis_other/ancestry/pca/merge_ref_samples.eigenvec",
        metadata="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam",
    output:
        plot="analysis_other/ancestry/pca/pca_plot.{population}.pdf"
    params:
        script="workflow/scripts/30_plot_pca.R"
    threads: 1
    log:
        "logs/ancestry/plot_pca/pca_plot_{population}.log"
    conda:
        "../env/r_ancestry.yml"
    shell:
        """
        Rscript \
            {params.script} \
            {input.matrix} \
            {output.plot} \
            NA \
            {input.metadata} \
            {wildcards.population} > {log} 2>&1
        """

# --------------- #
#    Admixture    #
# --------------- #

# Convert chromosome codes for ADMIXTURE (requires integer codes)
rule prepare_bed_for_admixture:
    input:
        pgen="analysis_other/ancestry/plink/merge_ref_samples.filtered.pgen",
        pvar="analysis_other/ancestry/plink/merge_ref_samples.filtered.pvar",
        psam="analysis_other/ancestry/plink/merge_ref_samples.filtered.psam"
    output:
        bed="analysis_other/ancestry/global/admixture/merged_cohort_numeric.bed",
        bim="analysis_other/ancestry/global/admixture/merged_cohort_numeric.bim",
        fam="analysis_other/ancestry/global/admixture/merged_cohort_numeric.fam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/prepare_bed_for_admixture.log"
    threads: 16
    shell:
        """
        mkdir -p analysis_other/ancestry/global/admixture
        
        # Convert chromosome codes to numeric format for ADMIXTURE
        plink2 \
            --pfile analysis_other/ancestry/plink/merge_ref_samples.filtered \
            --make-bed \
            --out analysis_other/ancestry/global/admixture/merged_cohort_numeric \
            --chr 1-22,X \
            --output-chr 26 \
            --threads {threads} \
            >{log} 2>&1
        """

rule run_admixture:
    input:
        bed="analysis_other/ancestry/global/admixture/merged_cohort_numeric.bed",
        bim="analysis_other/ancestry/global/admixture/merged_cohort_numeric.bim",
        fam="analysis_other/ancestry/global/admixture/merged_cohort_numeric.fam"
    output:
        touch("analysis_other/ancestry/global/admixture/results.done"),
        results="analysis_other/ancestry/global/admixture/merged_cohort_numeric.5.Q"
    conda:
        "../env/admixture.yml"
    log:
        "logs/ancestry/admixture.log"
    threads: 4
    shell:
        """
        mkdir -p analysis_other/ancestry/global/admixture
        
        # Store absolute path of log file before changing directory
        LOGFILE=$(realpath {log})
        
        # Run ADMIXTURE from the target directory (files are already in the right place)
        cd analysis_other/ancestry/global/admixture
        
        # Run ADMIXTURE for different K values
        for K in {{3..7}}; do
            echo "Running ADMIXTURE with K=$K" >>$LOGFILE 2>&1
            admixture --cv merged_cohort_numeric.bed $K -j{threads} >>$LOGFILE 2>&1
        done
        
        echo "ADMIXTURE analysis completed" >>$LOGFILE
        """

# --------------- #
#      iAdmix     #
# --------------- #

#############`###################################
# L O C A L  A N C E S T R Y   A N A L Y S I S #
################################################

# Prepeare Ref data for Local Ancestry
rule extract_ref_sample_list:
    input:
        vcf=config['1000G_phased_ref_vcf_biallelic']
    output:
        samples="analysis_other/ancestry/local/input/ref_sample_list.txt"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/extract_ref_sample_list.log"
    shell:
        """
        bcftools query -l {input.vcf} > {output.samples} 2>{log}
        """

rule create_sample_map:
    input:
        samples="analysis_other/ancestry/local/input/ref_sample_list.txt",
        panel=config['1000G_metadata_file']
    output:
        sample_map="analysis_other/ancestry/local/input/reference_sample_map.txt"
    params:
        pop_level="superpopulation"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/create_sample_map.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/local/input

        python3 workflow/scripts/25_create_ancestry_sample_map.py \
            --sample-list {input.samples} \
            --panel {input.panel} \
            --pop-level {params.pop_level} \
            --output {output.sample_map} \
            >{log} 2>&1
        """


# New rule to split rec VCF by chromosome for GnomiX
rule split_ref_vcf_by_chromosome:
    input:
        vcf=config['1000G_phased_ref_vcf_biallelic']
    output:
        vcf="analysis_other/ancestry/local/input/ref/1000G_phased.{chr}.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/split_vcf/ref_{chr}.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/local/input/ref
        
        bcftools view \
            --regions {wildcards.chr} \
            -Oz -o {output.vcf} \
            {input.vcf} \
            >{log} 2>&1
        
        tabix -p vcf {output.vcf} >>{log} 2>&1
        """

# Use the merged sample_ref VCF. Contains only intersection of SNPs in both QRY and REF
rule remove_ref_samples:
    input:
        vcf="analysis_other/ancestry/processed_vcf/merged_samples_ref.vcf.gz",
    output:
        vcf="analysis_other/ancestry/local/input/query/all_samples.ALL.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/remove_ref_samples.log"
    params:
        samples = lambda wc: ",".join(finished_samples)
    shell:
        """
        bcftools view \
            --samples "{params.samples}" \
            -Oz -o {output.vcf} \
            {input.vcf} \
            >{log} 2>&1
        tabix -p vcf {output.vcf} >>{log} 2>&1
        """
    
# Split VCF into chromosomes
rule split_vcf_by_chromosome:
    input:
        vcf="analysis_other/ancestry/local/input/query/all_samples.ALL.vcf.gz"
    output:
        vcf="analysis_other/ancestry/local/input/query/all_samples.{chr}.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/split_vcf/split_vcf_{chr}.log"
    shell:
        """      
        bcftools view \
            --regions {wildcards.chr} \
            -Oz -o {output.vcf} \
            {input.vcf} \
            >{log} 2>&1
        
        tabix -p vcf {output.vcf} >>{log} 2>&1
        """

# Modified rule for chromosome-specific genetic maps
rule extract_chromosome_genetic_map:
    input:
        gmap=config['1000G_genome_map']
    output:
        map="analysis_other/ancestry/local/genetic_map/{chr}.map"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/extract_chr_genetic_map/{chr}.log"
    run:
        with gzip.open(input.gmap, 'rt') as f_in, open(output.map, 'w') as f_out, open(log[0], 'w') as log_file:
            try:
                # Skip header
                header = f_in.readline()
                log_file.write(f"Processing chromosome {wildcards.chr}, skipping header: {header.strip()}\n")
                
                # Write data lines: position and cM (cumulative genetic distance)
                count = 0
                for line in f_in:
                    if line.strip():
                        fields = line.strip().split()
                        if len(fields) >= 3:
                            pos = int(fields[0])
                            cM = float(fields[2])  # Use cumulative genetic distance
                            f_out.write(f"{wildcards.chr}\t{pos}\t{cM}\n")
                            count += 1
                
                log_file.write(f"Successfully wrote {count} genetic map positions for chr{wildcards.chr}\n")
                
            except Exception as e:
                log_file.write(f"Error processing genetic map: {str(e)}\n")
                raise e


# Rule to create combined genetic map for all chromosomes
rule create_combined_genetic_map:
    input:
        gmaps=expand("analysis_other/ancestry/local/genetic_map/chr{chr}.map", chr=CHROMOSOMES)
    output:
        combined_map="analysis_other/ancestry/local/genetic_map/all_chr.map"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/rfmix/create_combined_genetic_map.log"
    shell:
        """
        # Concatenate and sort genetic map files lexicographically by chromosome and position
        cat {input.gmaps} | sort -k1,1V -k2,2n > {output.combined_map}
        echo "Combined genetic map created with $(wc -l < {output.combined_map}) positions" >>{log}
        """

# Rule to run RFMix by chromosome
rule run_rfmix_by_chromosome:
    input:
        query_vcf="analysis_other/ancestry/local/input/query/all_samples.{chr}.vcf.gz",
        reference_vcf="analysis_other/ancestry/local/input/ref/1000G_phased.{chr}.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/reference_sample_map.txt",
        genome_map="analysis_other/ancestry/local/genetic_map/{chr}.map",
    output:
        sis="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.sis.tsv",
        msp="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.msp.tsv",
        fb="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.fb.tsv",
        q="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.rfmix.Q"
    conda:
        "../env/rfmix.yml"
    resources:
        mem_gb=150
    log:
        "logs/ancestry/rfmix/{chr}/run_rfmix.log"
    threads: 6
    shell:
        """
        rfmix \
            -f {input.query_vcf} \
            -r {input.reference_vcf} \
            --chromosome={wildcards.chr} \
            --sample-map={input.sample_map} \
            --genetic-map={input.genome_map} \
            --output-basename=analysis_other/ancestry/local/rfmix/{wildcards.chr}/rfmix_{wildcards.chr} \
            --n-threads={threads} \
            --crf-spacing=0.001 \
            --rf-window-size=0.1 \
            >{log} 2>&1
        
        echo "RFMix whole genome analysis completed successfully" >>{log}
        """

rule combine_rfmix_results:
    input:
        msps=expand("analysis_other/ancestry/local/rfmix/chr{chr}/rfmix_chr{chr}.rfmix.Q", 
                   chr=[str(i) for i in range(1, 23)]) # chrX does not work (segfault, skipping for now)
    output:
        touch("analysis_other/ancestry/local/rfmix/results.done"),
    shell:
        """
        touch {output}
        """

# Rule to combine per-chromosome MSP files for each sample
rule combine_rfmix_msp_genome_wide:
    input:
        msp_files=expand("analysis_other/ancestry/local/rfmix/chr{chr}/rfmix_chr{chr}.msp.tsv",
                        chr=[str(i) for i in range(1, 23)])
    output:
        combined_msp="analysis_other/ancestry/local/rfmix/combined/rfmix_all_chr.msp.tsv"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/rfmix/combine_msp_genome_wide.log"
    shell:
        """
        python3 workflow/scripts/32_combine_rfmix_msp.py \
            --msp-files {input.msp_files} \
            --output {output.combined_msp} \
            >{log} 2>&1
        """

# Rule to convert RFMix MSP to Tagore BED format for each sample
rule rfmix_msp_to_tagore_bed:
    input:
        msp="analysis_other/ancestry/local/rfmix/combined/rfmix_all_chr.msp.tsv"
    output:
        bed="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.tagore.bed"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/rfmix/msp_to_bed_{sample}.log"
    shell:
        """
        python3 workflow/scripts/33_rfmix_to_tagore_bed.py \
            --msp {input.msp} \
            --sample {wildcards.sample} \
            --output {output.bed} \
            >{log} 2>&1
        """

# Rule to run Tagore for ancestry painting
rule run_tagore_ancestry_painting:
    input:
        bed="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.tagore.bed"
    output:
        plot="analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry.png"
    conda:
        "../env/tagore.yml"
    log:
        "logs/ancestry/rfmix/tagore_{sample}.log"
    params:
        prefix="analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry"
    threads: 1
    shell:
        """
        tagore \
            --input {input.bed} \
            --prefix {params.prefix} \
            --build hg38 \
            --oformat png \
            --force \
            --verbose \
            >{log} 2>&1
        """

# Rule to generate all Tagore plots
rule all_tagore_plots:
    input:
        expand("analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry.png",
               sample=finished_samples)
    output:
        touch("analysis_other/ancestry/local/rfmix/tagore_plots.done")



# Modified rule to train GnomiX model per chromosome
rule gnomix_training_by_chr:
#    <genetic_map_file> is the genetic map file. It's a .tsv file with 3 columns; chromosome number, SNP physical position and SNP genetic position. There should be no headers unless they start with "#". See example in the demo/data/ folder.
#    <reference_file> is a .vcf or .vcf.gz file containing the reference haplotypes (in any order)
#    <sample_map_file> is a sample map file matching reference samples to their respective reference populations
    input:
        query_vcf="analysis_other/ancestry/local/input/query/all_samples.{chr}.vcf.gz",
        reference_vcf="analysis_other/ancestry/local/input/ref/1000G_phased.{chr}.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/reference_sample_map.txt",
        genome_map="analysis_other/ancestry/local/genetic_map/{chr}.map",
    output:
        model="analysis_other/ancestry/local/gnomix/model/{chr}/models/model_chm_{chr}/model_chm_{chr}.pkl",
        msp="analysis_other/ancestry/local/gnomix/model/{chr}/query_results.msp"
    conda:
        "../env/gnomix.yml"
    log:
        "logs/ancestry/gnomix/run_gnomix_with_training_{chr}.log"
    resources:
        mem_gb=200
    params:
        gnomix=config['gnomix'],
        config=config['gnomix_config'],
    threads: 40
    # python3 gnomix.py <query_file> <output_folder> <chr_nr> <phase> <genetic_map_file> <reference_file> <sample_map_file>
    shell:
        """
        {params.gnomix} \
            {input.query_vcf} \
            $(dirname {output.msp}) \
            {wildcards.chr} \
            True \
            {input.genome_map} \
            {input.reference_vcf} \
            {input.sample_map} \
            {params.config} > {log} 2>&1
        """

rule combine_gnomix_results:
    input:
        msps=expand("analysis_other/ancestry/local/gnomix/model/{chr}/query_results.msp",
                   chr=["chr" + str(i) for i in range(1, 23)] + ["chrX"])
    output:
        touch("analysis_other/ancestry/local/gnomix/results.done"),
        combined="analysis_other/ancestry/local/gnomix/results/gnomix_lai.msp"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/combine_gnomix_results.log"
    run:
        log_file = open(log[0], "w")
        try:
            dfs = []
            for f in input.msps:
                try:
                    with open(f) as fh:
                        header_line = fh.readline()
                    df = pd.read_csv(f, sep='\t', skiprows=1)
                    dfs.append(df)
                    log_file.write(f"Read {f}: {len(df)} rows\n")
                except Exception as e:
                    log_file.write(f"Error reading {f}: {str(e)}\n")
            
            if dfs:
                combined = pd.concat(dfs, ignore_index=True)
                combined = combined.sort_values(['chm', 'spos', 'epos'])
                os.makedirs(os.path.dirname(output.combined), exist_ok=True)
                combined.to_csv(output.combined, sep='\t', index=False)
                log_file.write(f"Combined {len(dfs)} files, {len(combined)} total rows\n")
            else:
                with open(output.combined, 'w') as f:
                    f.write("# No valid GnomiX results\n")
        except Exception as e:
            log_file.write(f"Error: {str(e)}\n")
            raise e
        finally:
            log_file.close()

