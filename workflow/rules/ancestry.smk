# Global variables for ancestry analysis
ANCESTRY_TOOLS = ["iadmix"]

LOCAL_ANCESTRY_TOOLS = ["gnomix"]
#LOCAL_ANCESTRY_TOOLS = ["rfmix"]

THOUSAND_G_POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

rule all_ancestry:
    input:
        # Merged and QC'd variant set
        "analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        "analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        "analysis_other/ancestry/plink/merged/merged_cohort.psam",
        # 1000G reference data (filtered)
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.pgen",
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.pvar",
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam",
        # Sample data (filtered)
        "analysis_other/ancestry/plink/samples/samples_filtered.pgen",
        "analysis_other/ancestry/plink/samples/samples_filtered.pvar",
        "analysis_other/ancestry/plink/samples/samples_filtered.psam",
        # Global ancestry results
        expand("analysis_other/ancestry/global/iadmix/freq/freq_{population}.frq", population=THOUSAND_G_POPS),
        expand("analysis_other/ancestry/global/{tool}/results.done", tool=ANCESTRY_TOOLS),
        # Local ancestry results  
        expand("analysis_other/ancestry/local/{tool}/results.done", tool=LOCAL_ANCESTRY_TOOLS),
        # PCA plots for each population
        expand("analysis_other/ancestry/pca/pca_plot.{population}.pdf", population=["POP", "SUPERPOP"]),
        # RFMIX tagore plots
        #"analysis_other/ancestry/local/rfmix/tagore_plots.done"
        # Summary reports
        #"analysis_other/ancestry/reports/ancestry_summary.html"



# Part A: Preprocess 1000G ref data

# Step A1: Add chr prefix and basic processing
rule add_chr_prefix_1000g:
    output:
        vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_prefixed.vcf.gz"),
        tbi=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_prefixed.vcf.gz.tbi")
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_add_chr_prefix/{chr}.log"
    shell:
        """
        # Add chr prefix using bcftools annotate
        bcftools annotate \
            --rename-chrs <(echo "{wildcards.chr} chr{wildcards.chr}") \
            --output-type z \
            --output {output.vcf} \
            data/ref/variant_sets/1000G/ALL.chr{wildcards.chr}.*.vcf.gz \
            2>{log}
        
        tabix -p vcf {output.vcf} 2>>{log}
        """

# Step A2: Normalize and split multiallelic variants
rule normalize_and_split_multiallelic:
    input:
        vcf="analysis_other/ancestry/processed_vcf/1000G/chr{chr}_prefixed.vcf.gz",
    output:
        vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_normalized.vcf.gz"),
        tbi=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_normalized.vcf.gz.tbi")
    params:
        ref=config['ref_grch37']
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_normalize_chr/{chr}.log"
    shell:
        """
        # Normalize and split multiallelic variants
        bcftools norm \
            --multiallelics -both \
            --fasta-ref {params.ref} \
            --check-ref w \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} \
            2>{log}
        
        tabix -p vcf {output.vcf} 2>>{log}
        """

# Step A3: Filter out problematic variants before liftover
rule filter_for_liftover:
    input:
        vcf="analysis_other/ancestry/processed_vcf/1000G/chr{chr}_normalized.vcf.gz"
    output:
        vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_filtered.vcf.gz"),
        tbi=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_filtered.vcf.gz.tbi")
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_filter_chr/{chr}.log"
    shell:
        """
        # Filter out structural variants and keep only SNPs and simple indels, remove half-calls
        bcftools view \
            --types snps,indels \
            --exclude 'ALT~"<" || REF~"<" || ALT~">" || REF~">" || GT="mis"' \
            --min-alleles 2 \
            --max-alleles 2 \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} \
            >{log} 2>&1
        
        tabix -p vcf {output.vcf} >>{log} 2>&1
        """

# Step A4 Liftover coordinates with CrossMap (produces unsorted VCF)
rule crossmap_liftover:
    input:
        vcf="analysis_other/ancestry/processed_vcf/1000G/chr{chr}_filtered.vcf.gz"
    output:
        vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.unsorted.vcf")
    params:
        chain=config['liftover_chain_GRCh37_to_T2T'],
        ref=config['ref']
    conda:
        "../env/liftover.yml"
    log:
        "logs/ancestry/ref_liftover_T2T_chr/{chr}.log"
    shell:
        """
        CrossMap.py vcf {params.chain} \
            {input.vcf} \
            {params.ref} \
            {output.vcf} >>{log} 2>&1
        """

# Step A5: Fix VCF header
rule add_contigs_to_header:
    input:
        vcf="analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.unsorted.vcf",
        fai=config['ref_fai']
    output:
        vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.headerfix.vcf")
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_add_contigs_chr/{chr}.log"
    shell:
        """
        bcftools reheader --fai {input.fai} -o {output.vcf} {input.vcf} 2>{log}
        """

# Step A6: Sort, compress, and index the lifted VCF
rule sort_compress_index_lifted_vcf:
    input:
        vcf="analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.headerfix.vcf"
    output:
        lifted_vcf=temp("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.vcf.gz")
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_sort_compress_index_chr/{chr}.log"
    shell:
        """
        bcftools sort {input.vcf} -Oz -o {output.lifted_vcf} 2>{log}
        tabix -p vcf {output.lifted_vcf} 2>>{log}
        """

# Step A7: Concatenate 1000G chromosomes into a merged VCF
rule concat_1000g_vcfs:
    input:
        vcfs=expand("analysis_other/ancestry/processed_vcf/1000G/chr{chr}_T2T.vcf.gz", chr=CHROMOSOMES)
    output:
        merged_vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/ref_concat_1000g_vcfs.log"
    shell:
        """
        bcftools concat \
            {input.vcfs} \
            --output-type z \
            --output {output.merged_vcf} \
            >{log} 2>&1
        tabix -p vcf {output.merged_vcf}
        """

# Step A8: Convert metadata to PSAM, using panel file for population info
rule vcf_1000G_tsv_to_psam:
    input:
        metadata="data/ref/variant_sets/1000G/integrated_call_samples_v3.20200731.ALL.ped",
        panel="data/ref/variant_sets/1000G/integrated_call_samples_v3.20130502.ALL.panel",
        vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz"
    output:
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.psam"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/ref_metadata_to_psam.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/plink/reference

        python3 workflow/scripts/26_convert_metadata_to_psam.py \
            --vcf {input.vcf} \
            --metadata {input.metadata} \
            --panel {input.panel} \
            --psam {output.psam} \
            >{log} 2>&1
        """

# Step A9: Convert VCF to PGEN format
rule vcf_to_plink_1000g:
    input:
        merged_vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz",
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.psam"
    output:
        pgen="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.pgen",
        pvar="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.pvar",
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/vcf_to_plink_1000g.log"
    threads: 32
    shell:
        """
        plink2 \
            --vcf {input.merged_vcf} \
            --make-pgen \
            --psam {input.psam} \
            --out analysis_other/ancestry/plink/reference/1000G_phase3_T2T \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --threads {threads} \
            >{log} 2>&1
        """


# Process samples data

# Step B1: Merge sample VCFs 
rule merge_sample_vcfs:
    input:
        vcfs=lambda wc: expand("assembly/variants/{sample}/phased_verkko/small_variants.vcf.gz", 
                              sample=[s for s in finished_samples if s in asm_samples])
    output:
        merged="analysis_other/ancestry/processed_vcf/samples/merged_samples.vcf.gz",
        list="analysis_other/ancestry/processed_vcf/samples/vcf_list.txt"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/sample_merge_vcfs.log"
    shell:
        """
        # Create list of VCF files
        echo {input.vcfs} | tr ' ' '\n' > {output.list}
        
        # Normalize and merge VCFs
        bcftools merge \
            --merge none \
            --file-list {output.list} \
            --output-type z \
            --output {output.merged} \
            >{log} 2>&1
        
        tabix -p vcf {output.merged}
        """

# Step B2: filter out problematic variants from sample VCFs
rule filter_sample_vcf:
    input:
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples.vcf.gz"
    output:
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/sample_filter_vcf.log"
    shell:
        """
        # First normalize and split multiallelic sites, then filter out large variants, structural variants
        bcftools norm \
            --multiallelics -both \
            --output-type u \
            {input.vcf} 2>{log} | \
        bcftools view \
            --types snps,indels \
            --exclude 'ALT~"\\*" || (strlen(REF)>50 || strlen(ALT)>50)' \
            --min-alleles 2 \
            --max-alleles 2 \
            --output-type z \
            --output {output.vcf} \
            >>{log} 2>&1
        
        tabix -p vcf {output.vcf} 2>>{log}
        """

# Step B3: Convert sample VCF to PGEN format
rule samples_to_plink:
    input:
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz",
        metadata="data/samples_pedigree.tsv"
    output:
        pgen="analysis_other/ancestry/plink/samples/samples.pgen",
        pvar="analysis_other/ancestry/plink/samples/samples.pvar",
        psam="analysis_other/ancestry/plink/samples/samples.psam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/samples_to_plink.log"
    threads: 32
    shell:
        """
        plink2 \
            --vcf {input.vcf} \
            --make-pgen \
            --out analysis_other/ancestry/plink/samples/samples \
            --psam {input.metadata} \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --set-missing-var-ids @:# \
            --threads {threads} \
            >{log} 2>&1
        """


# Step C1: Filter datasets with QC filters using PGEN format
rule qc_filter_plink:
    input:
        pgen="analysis_other/ancestry/plink/{dataset}/{prefix}.pgen",
        pvar="analysis_other/ancestry/plink/{dataset}/{prefix}.pvar",
        psam="analysis_other/ancestry/plink/{dataset}/{prefix}.psam"
    output:
        pgen="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.pgen",
        pvar="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.pvar",
        psam="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.psam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/qc_filter_{dataset}_{prefix}.log"
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

        echo "Applying QC filters to {wildcards.dataset}/{wildcards.prefix}..." >{log}
        
        plink2 \
            --pfile analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix} \
            --maf $MAF_THRESHOLD \
            --geno $GENO_THRESHOLD \
            --hwe $HWE_THRESHOLD \
            --make-pgen \
            --out analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_filtered \
            --threads {threads} \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 1000 \
            --sort-vars \
            >{log} 2>&1
        
        # Report filtering statistics
        echo "Filtering completed for {wildcards.dataset}/{wildcards.prefix}" >>{log}
        echo "Variants after filtering: $(wc -l < analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_filtered.pvar | tail -n +2)" >>{log}
        echo "Samples after filtering: $(wc -l < analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_filtered.psam | tail -n +2)" >>{log}
        """

# Step C2: Convert PGEN to BED format for merging (works for both reference and samples)
rule convert_pgen_to_bed_for_merge:
    input:
        pgen="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.pgen",
        pvar="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.pvar",
        psam="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.psam"
    output:
        bed=temp("analysis_other/ancestry/plink/{dataset}/{prefix}_old.bed"),
        bim=temp("analysis_other/ancestry/plink/{dataset}/{prefix}_old.bim"),
        fam=temp("analysis_other/ancestry/plink/{dataset}/{prefix}_old.fam")
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/convert_to_bed_{dataset}_{prefix}.log"
    threads: 16
    shell:
        """
        mkdir -p analysis_other/ancestry/plink/merged
        
        # Convert PGEN to BED format using PLINK2
        plink2 \
            --pfile analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_filtered \
            --make-bed \
            --out analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_old \
            --threads {threads} \
            >{log} 2>&1
        """

# Step C3: Merge BED files using PLINK 1.9
rule merge_bed_files_plink:
    input:
        ref_bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.bed",
        ref_bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.bim", 
        ref_fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.fam",
        sample_bed="analysis_other/ancestry/plink/samples/samples_old.bed",
        sample_bim="analysis_other/ancestry/plink/samples/samples_old.bim",
        sample_fam="analysis_other/ancestry/plink/samples/samples_old.fam"
    output:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam",
        merge_report="analysis_other/ancestry/qc/merge_report.txt"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/merge_bed_files.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/qc
        
        # Create merge list for PLINK 1.9
        echo "analysis_other/ancestry/plink/samples/samples_old" > analysis_other/ancestry/plink/merged/merge_list.txt
        
        # Merge using PLINK 1.9 with --merge
        plink \
            --bfile analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old \
            --merge-list analysis_other/ancestry/plink/merged/merge_list.txt \
            --out analysis_other/ancestry/plink/merged/merged_cohort \
            --allow-extra-chr \
            --threads {threads} \
            >{log} 2>&1
        
        # Generate merge report
        echo "Merge completed on $(date)" > {output.merge_report}
        
        # Count samples and variants from BED/BIM/FAM files
        if [[ -f {output.fam} ]]; then
            echo "Final sample count: $(wc -l < {output.fam})" >> {output.merge_report}
        fi
        
        if [[ -f {output.bim} ]]; then
            echo "Final variant count: $(wc -l < {output.bim})" >> {output.merge_report}
        fi
        
        # Check for merge conflicts
        if [[ -f analysis_other/ancestry/plink/merged/merged_cohort_temp.missnp ]]; then
            echo "Variants with merge conflicts: $(wc -l < analysis_other/ancestry/plink/merged/merged_cohort_temp.missnp)" >> {output.merge_report}
        fi
        
        # Clean up merge list
        rm -f analysis_other/ancestry/plink/merged/merge_list.txt
        """

# Step C4: Convert merged BED back to PGEN format
rule convert_bed_to_pgen_final:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim", 
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam",
        merge_report="analysis_other/ancestry/qc/merge_report.txt"
    output:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort_no_pop.psam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/convert_to_pgen_final.log"
    threads: 32
    shell:
        """
        # Convert merged BED back to PGEN format using PLINK2
        plink2 \
            --bfile analysis_other/ancestry/plink/merged/merged_cohort \
            --make-pgen \
            --out analysis_other/ancestry/plink/merged/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        
        # Move PSAM to indicate no population info yet
        mv analysis_other/ancestry/plink/merged/merged_cohort.psam {output.psam}

        # Update merge report with final PGEN file statistics
        echo "" >> {input.merge_report}
        echo "Conversion to PGEN completed on $(date)" >> {input.merge_report}
        echo "Final PGEN sample count: $(tail -n +2 {output.psam} | wc -l)" >> {input.merge_report}
        echo "Final PGEN variant count: $(tail -n +2 {output.pvar} | wc -l)" >> {input.merge_report}
        """

# Step C5: Restore population information to merged PSAM file
rule restore_population_info_to_merged_psam:
    input:
        merged_psam="analysis_other/ancestry/plink/merged/merged_cohort_no_pop.psam",
        ref_psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam",
        sample_psam="analysis_other/ancestry/plink/samples/samples_filtered.psam"
    output:
        enhanced_psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/restore_population_info.log"
    shell:
        """
        python3 workflow/scripts/27_restore_population_info.py \
            --merged-psam {input.merged_psam} \
            --ref-psam {input.ref_psam} \
            --sample-psam {input.sample_psam} \
            --output {output.enhanced_psam} \
            >{log} 2>&1
        """

# Part D: Ancestry analysis

# Step D1: PCA for population structure using PGEN format
rule run_pca:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar", 
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        eigenval="analysis_other/ancestry/pca/merged_cohort.eigenval",
        eigenvec="analysis_other/ancestry/pca/merged_cohort.eigenvec"
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
            --pfile analysis_other/ancestry/plink/merged/merged_cohort \
            --indep-pairwise 50 10 0.2 \
            --out analysis_other/ancestry/pca/ld_pruned \
            --threads {threads} \
            >{log} 2>&1
        
        # Run PCA
        plink2 \
            --pfile analysis_other/ancestry/plink/merged/merged_cohort \
            --extract analysis_other/ancestry/pca/ld_pruned.prune.in \
            --pca 20 \
            --out analysis_other/ancestry/pca/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        """


rule plot_pca:
    input:
        matrix="analysis_other/ancestry/pca/merged_cohort.eigenvec",
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


# Step D2a-1: Create population-specific keep files
rule create_population_keep_file:
    input:
        ref_psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam"
    output:
        keep_file="analysis_other/ancestry/global/iadmix/keep_{population}.txt"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/iadmix/create_keep_file_{population}.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/global/iadmix
        
        POP={wildcards.population}
        echo "Creating keep file for population: $POP" >{log}
        
        # Create population-specific keep file from PSAM (column 8 is SUPERPOP, skip header with #)
        awk -v pop="$POP" '!/^#/ && $8 == pop {{print $1, $2}}' {input.ref_psam} > {output.keep_file}
        
        # Debug: Show first few lines of PSAM file and what we're looking for
        echo "Looking for population: $POP" >>{log}
        echo "PSAM header:" >>{log}
        head -1 {input.ref_psam} >>{log}
        echo "Sample PSAM lines:" >>{log}
        head -5 {input.ref_psam} | tail -4 >>{log}
        echo "Generated keep file content:" >>{log}
        head -5 {output.keep_file} >>{log}
        
        # Report count of samples in keep file
        SAMPLE_COUNT=$(wc -l < {output.keep_file})
        echo "Found $SAMPLE_COUNT samples for population $POP" >>{log}
        """

# Step D2a-2: Calculate allele frequencies using population-specific keep files
rule calculate_plink_frequencies:
    input:
        ref_bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.bed",
        ref_bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.bim", 
        ref_fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old.fam",
        keep_file="analysis_other/ancestry/global/iadmix/keep_{population}.txt"
    output:
        freq_file="analysis_other/ancestry/global/iadmix/freq/freq_{population}.frq"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/iadmix/calculate_freq_{population}.log"  
    threads: 4  
    shell:
        """
        mkdir -p analysis_other/ancestry/global/iadmix/freq
        
        POP={wildcards.population}
        echo "Calculating frequencies for population: $POP" >{log}
        
        # Check if keep file has any samples
        SAMPLE_COUNT=$(wc -l < {input.keep_file})
        echo "Using $SAMPLE_COUNT samples for population $POP" >>{log}
        
        if [[ $SAMPLE_COUNT -gt 0 ]]; then
            # Calculate frequencies using PLINK 1.9
            plink \
                --bfile analysis_other/ancestry/plink/reference/1000G_phase3_T2T_old \
                --keep {input.keep_file} \
                --freq \
                --out analysis_other/ancestry/global/iadmix/freq_{wildcards.population} \
                --allow-extra-chr \
                >>{log} 2>&1
            
            # Check if freq file exists and has correct name
            if [[ -f analysis_other/ancestry/global/iadmix/freq_{wildcards.population}.frq ]]; then
                cp analysis_other/ancestry/global/iadmix/freq_{wildcards.population}.frq {output.freq_file}
            else
                echo "PLINK frequency calculation failed for $POP" >>{log}
                touch {output.freq_file}
            fi
        else
            echo "No samples found for population $POP, creating empty frequency file" >>{log}
            touch {output.freq_file}
        fi
        """

# Step D2b: Convert PLINK frequency output to iAdmix format
rule convert_to_iadmix_freq_format:
    input:
        plink_freq=expand("analysis_other/ancestry/global/iadmix/freq/freq_{population}.frq", population=THOUSAND_G_POPS),  # Fixed extension
        pvar="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.pvar",
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam"
    output:
        freq_file="analysis_other/ancestry/global/iadmix/freq/reference_frequencies.txt"
    conda:
        "../env/py_report.yml"  # Use Python environment
    log:
        "logs/ancestry/iadmix/convert_iadmix_freq.log"
    shell:
        """
        # Convert PLINK frequency output to iAdmix format
        python3 workflow/scripts/28_create_iadmix_freq.py \
            --plink-freq {input.plink_freq} \
            --pvar {input.pvar} \
            --psam {input.psam} \
            --output {output.freq_file} \
            >{log} 2>&1
        """

# Function to get sample IDs from PSAM file
def get_query_sample_ids():
    import pandas as pd
    psam_file = "analysis_other/ancestry/plink/merged/merged_cohort.psam"
    try:
        psam_df = pd.read_csv(psam_file, sep='\t')
        psam_df.columns = [col.lstrip('#') for col in psam_df.columns]
        # Get only query samples (not 1000G reference)
        if 'SUPERPOP' in psam_df.columns:
            query_samples = psam_df[psam_df['SUPERPOP'] == 'QUERY']['IID'].tolist()
        else:
            # If no SUPERPOP column, get all samples (fallback)
            query_samples = psam_df['IID'].tolist()
        return query_samples
    except:
        # Fallback if file doesn't exist yet
        return []

# Get sample IDs (this will be populated after the merged PSAM file is created)
QUERY_SAMPLES = get_query_sample_ids()

# Step D2c: Export sample genotypes for individual processing
rule export_iadmix_plink_genotypes_individual:
    input:
        pgen="analysis_other/ancestry/plink/samples/samples_filtered.pgen",
        pvar="analysis_other/ancestry/plink/samples/samples_filtered.pvar",
        psam="analysis_other/ancestry/plink/samples/samples_filtered.psam"
    output:
        raw_file="analysis_other/ancestry/global/iadmix/geno/samples_raw.traw"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/iadmix/export_plink_geno_individual.log"
    threads: 16
    shell:
        """
        mkdir -p analysis_other/ancestry/global/iadmix/geno
        
        # Export genotypes in A-transpose format (samples as columns, variants as rows)
        plink2 \
            --pfile analysis_other/ancestry/plink/samples/samples_filtered \
            --export A-transpose \
            --out analysis_other/ancestry/global/iadmix/geno/samples_raw \
            --threads {threads} \
            >{log} 2>&1
        """

# Step D2d: Convert PLINK genotypes to iAdmix format for individual samples
rule convert_to_iadmix_geno_format_individual:
    input:
        plink_raw="analysis_other/ancestry/global/iadmix/geno/samples_raw.traw"
    output:
        geno_file="analysis_other/ancestry/global/iadmix/geno/sample_{sample}_genotypes.txt"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/iadmix/convert_geno_{sample}.log"
    shell:
        """
        # Convert to iAdmix genotype format for specific sample
        python3 workflow/scripts/29_create_iadmix_geno.py \
            --plink-raw {input.plink_raw} \
            --sample {wildcards.sample} \
            --output {output.geno_file} \
            >{log} 2>&1
        """

# Step D2e: Run iAdmix for individual samples
rule run_iadmix_individual:
    input:
        freq_file="analysis_other/ancestry/global/iadmix/freq/reference_frequencies.txt",
        geno_file="analysis_other/ancestry/global/iadmix/geno/sample_{sample}_genotypes.txt"
    output:
        results="analysis_other/ancestry/global/iadmix/results/sample_{sample}_ancestry.txt"
    log:
        "logs/ancestry/iadmix/run_{sample}.log"
    threads: 4
    params:
        wd = workflow.basedir
    shell:  
        """
        docker run --rm \
            -u $(id -u):$(id -g) \
            -v $(realpath .):/workdir \
            -w /workdir \
            caspargross/iadmix \
            python /usr/src/app/runancestry.py \
            --freq {input.freq_file} \
            --geno {input.geno_file} \
            --out analysis_other/ancestry/global/iadmix/results/sample_{wildcards.sample} \
            --cores {threads} \
            --path /usr/src/app/ \
            >{log} 2>&1
        
        rm {input.geno_file}.ancestry.input   # Clean up genotype file after processing

        if [[ -f analysis_other/ancestry/global/iadmix/results/sample_{wildcards.sample} ]]; then
            # Extract ancestry proportions and create formatted output
            echo -e "{wildcards.sample}\t$(tail -1 analysis_other/ancestry/global/iadmix/results/sample_{wildcards.sample})" > {output.results}
        else
            # Create empty result if analysis failed
            echo "iAdmix analysis failed for {wildcards.sample}" >>{log}
        fi
        """

# Step D2f: Combine all individual iAdmix results
rule combine_iadmix_individual_results:
    input:
        lambda wildcards: expand("analysis_other/ancestry/global/iadmix/results/sample_{sample}_ancestry.txt", 
                                sample=QUERY_SAMPLES) if QUERY_SAMPLES else []
    output:
        combined="analysis_other/ancestry/global/iadmix/all_samples_ancestry_individual.txt",
        touch_file=touch("analysis_other/ancestry/global/iadmix/results.done")
    log:
        "logs/ancestry/iadmix/combine_individual_results.log"
    run:
        import os
        
        # Create output directory
        os.makedirs("analysis_other/ancestry/global/iadmix", exist_ok=True)
        
        # Create header
        with open(output.combined, 'w') as f:
            f.write("Sample_ID\tAFR\tAMR\tEAS\tEUR\tSAS\n")
        
        # Combine results
        if input:
            shell("cat {input} >> {output.combined} 2>{log}")
            with open(str(log), 'a') as f:
                f.write(f"Combined ancestry results for {len(input)} sample files\\n")
        else:
            with open(output.combined, 'a') as f:
                f.write("No samples processed\\n")
            with open(str(log), 'a') as f:
                f.write("No sample results to combine\\n")
        
        with open(str(log), 'a') as f:
            f.write("Individual iAdmix analysis completed\n")

# Update the main iadmix rule to use individual results
rule run_iadmix:
    input:
        individual_results="analysis_other/ancestry/global/iadmix/all_samples_ancestry_individual.txt"
    output:
        results="analysis_other/ancestry/global/iadmix/results/admixture_proportions.txt"
    log:
        "logs/ancestry/iadmix_final.log"
    shell:
        """
        # Copy individual results to main output location
        cp {input.individual_results} {output.results} 2>{log}
        echo "iAdmix individual sample analysis completed" >>{log}
        """


# Step D2a: Convert chromosome codes for ADMIXTURE (requires integer codes)
rule prepare_bed_for_admixture:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam"
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
            --bfile analysis_other/ancestry/plink/merged/merged_cohort \
            --make-bed \
            --out analysis_other/ancestry/global/admixture/merged_cohort_numeric \
            --chr 1-22,X \
            --output-chr 26 \
            --threads {threads} \
            >{log} 2>&1
        """

# Updated ADMIXTURE rule to use numeric chromosome codes
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

# Prepare data for local ancestry analysis - Step 1: Convert PGEN to VCF
rule convert_pgen_to_vcf:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        phased_vcf="analysis_other/ancestry/local/input/merged_cohort_phased.vcf.gz"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/convert_pgen_to_vcf.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/local/input
        
        # Convert PGEN to VCF format
        plink2 \
            --pfile analysis_other/ancestry/plink/merged/merged_cohort \
            --export vcf bgz \
            --out analysis_other/ancestry/local/input/merged_cohort_phased \
            --threads {threads} \
            >{log} 2>&1
        """

# Prepare data for local ancestry analysis - Step 2: Create sample mapping
rule create_sample_map:
    input:
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam"
    output:
        sample_map="analysis_other/ancestry/local/input/sample_map.txt"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/create_sample_map.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/local/input
        
        # Create sample mapping file for RFMix
        python3 workflow/scripts/25_create_ancestry_sample_map.py \
            --psam {input.psam} \
            --output {output.sample_map} \
            >{log} 2>&1
        """

# New rule to split VCF by chromosome for GnomiX
rule split_vcf_by_chromosome_gnomix:
    input:
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz"
    output:
        vcf="analysis_other/ancestry/local/input/samples.{chr}.vcf.gz"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/gnomix/split_vcf/{chr}.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/local/gnomix/input
        
        bcftools view \
            --regions {wildcards.chr} \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} \
            >{log} 2>&1
            
        bcftools index --tbi {output.vcf} 2>>{log}
        """

# Modified rule for chromosome-specific genetic maps
rule extract_chromosome_genetic_map:
    input:
        gmap="data/ref/phasing_T2T/resources/recombination_maps/t2t_native_scaled_maps/{chr}.t2t.scaled.gmap.gz"
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

# RFMix local ancestry analysis rules

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
rule run_rfmix_whole_genome:
    input:
        query_vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz",
        reference_vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt",
        genetic_map="analysis_other/ancestry/local/genetic_map/all_chr.map"
    output:
        sis="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.sis.tsv",
        msp="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.msp.tsv",
        fb="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.fb.tsv",
        q="analysis_other/ancestry/local/rfmix/{chr}/rfmix_{chr}.rfmix.Q"
    conda:
        "../env/rfmix.yml"
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
            --genetic-map={input.genetic_map} \
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

# Rule to convert RFMix MSP to BED format for each sample
rule rfmix_msp_to_bed:
    input:
        msp="analysis_other/ancestry/local/rfmix/combined/rfmix_all_chr.msp.tsv"
    output:
        hap0="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.hap0.bed",
        hap1="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.hap1.bed",
        config="analysis_other/ancestry/local/rfmix/tagore_input/{sample}_tagore.conf"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/rfmix/msp_to_bed_{sample}.log"
    shell:
        """
        python3 workflow/scripts/33_rfmix_to_tagore_bed.py \
            --msp {input.msp} \
            --sample {wildcards.sample} \
            --output-prefix analysis_other/ancestry/local/rfmix/tagore_input/{wildcards.sample} \
            --create-config \
            >{log} 2>&1
        """

# Rule to create chromosome sizes file for Tagore
rule create_chrom_sizes:
    input:
        fai=config['ref_fai']
    output:
        sizes="analysis_other/ancestry/local/rfmix/tagore_input/chrom.sizes"
    shell:
        """
        # Extract chromosome sizes for autosomes only
        awk '$1 ~ /^chr[0-9]+$/ {{print $1"\\t"$2}}' {input.fai} | \
            sort -V > {output.sizes}
        """

# Rule to run Tagore for ancestry painting
rule run_tagore_ancestry_painting:
    input:
        hap0="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.hap0.bed",
        hap1="analysis_other/ancestry/local/rfmix/tagore_input/{sample}.hap1.bed",
        config="analysis_other/ancestry/local/rfmix/tagore_input/{sample}_tagore.conf",
        sizes="analysis_other/ancestry/local/rfmix/tagore_input/chrom.sizes"
    output:
        plot="analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry.png"
    conda:
        "../env/tagore.yml"
    log:
        "logs/ancestry/rfmix/tagore_{sample}.log"
    params:
        output_prefix="analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry"
    threads: 1
    shell:
        """
        tagore \
            --input {input.hap0} {input.hap1} \
            --genome {input.sizes} \
            --color {input.config} \
            --output {params.output_prefix} \
            --prefix {wildcards.sample} \
            --height 800 \
            --width 1200 \
            --format png \
            >{log} 2>&1
        """

# Rule to generate all Tagore plots
rule all_tagore_plots:
    input:
        expand("analysis_other/ancestry/local/rfmix/plots/{sample}_ancestry.png",
               sample=QUERY_SAMPLES)
    output:
        touch("analysis_other/ancestry/local/rfmix/tagore_plots.done")

# Update the all_ancestry rule to include Tagore plots
# Modify the existing rule by adding to the input list:
# "analysis_other/ancestry/local/rfmix/tagore_plots.done"


# Modified rule to train GnomiX model per chromosome
rule run_gnomix_with_training_by_chr:
    input:
        query_vcf=ancient("analysis_other/ancestry/local/input/samples.{chr}.vcf.gz"),
        reference_vcf="analysis_other/ancestry/processed_vcf/1000G/{chr}_T2T.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt",
        genetic_map="analysis_other/ancestry/local/genetic_map/{chr}.map"
    output:
        model="analysis_other/ancestry/local/gnomix/{chr}/models/1000G_T2T_chm_{chr}/1000G_T2T_chm_{chr}.pkl",
        msp="analysis_other/ancestry/local/gnomix/{chr}/query_results.msp"
    conda:
        "../env/gnomix.yml"
    log:
        "logs/ancestry/gnomix/run_gnomix_with_training_{chr}.log"
    params:
        gnomix=config['gnomix'],
        config=config['gnomix_config'],
    threads: 1
    # python3 gnomix.py <query_file> <output_folder> <chr_nr> <phase> <genetic_map_file> <reference_file> <sample_map_file>
    shell:
        """
        
        # Use script command to capture all output including subprocess calls
        {params.gnomix} \
            {input.query_vcf} \
            $(dirname {output.msp}) \
            {wildcards.chr} \
            True \
            {input.genetic_map} \
            {input.reference_vcf} \
            {input.sample_map} \
            {params.config} > {log} 2>&1
        """

# Modified rule to train GnomiX model per chromosome
rule run_gnomix_by_chr:
    input:
        model=ancient("analysis_other/ancestry/local/gnomix/{chr}/models/1000G_T2T_chm_{chr}/1000G_T2T_chm_{chr}.pkl"),
        query_vcf="analysis_other/ancestry/local/input/samples.{chr}.vcf.gz",
    output:
        msp="analysis_other/ancestry/local/gnomix/{chr}/query_results.msp"
    conda:
        "../env/gnomix.yml"
    log:
        "logs/ancestry/gnomix/run_gnomix_with_training_{chr}.log"
    params:
        gnomix=config['gnomix'],
        config=config['gnomix_config'],
    threads: 1
        #$ python3 gnomix.py <query_file> <output_folder> <chr_nr> <phase> <path_to_model>     
    shell:
        """
        mkdir -p analysis_other/ancestry/local/gnomix/trained_model
        
        # Train GnomiX model for this chromosome
        {params.gnomix} \
            {input.query_vcf} \
            $(dirname {output.msp}) \
            {wildcards.chr} \
            True \
            $(dirname {input.model} \
            >{log} 2>&1
        """

ruleorder:  run_gnomix_with_training_by_chr > run_gnomix_by_chr

# Aggregate rule to combine results from all chromosomes
rule combine_gnomix_results:
    input:
        msps=expand("analysis_other/ancestry/local/gnomix/chr{chr}/query_results.msp", 
                   chr=[str(i) for i in range(1, 23)] + ["X"])
    output:
        touch("analysis_other/ancestry/local/gnomix/results.done"),
        combined="analysis_other/ancestry/local/gnomix/results/gnomix_lai.msp"
    conda:
        "../env/py_report.yml"
    log:
        "logs/ancestry/combine_gnomix_results.log"
    run:
        # Set up logging
        log_file = open(log[0], "w")
        
        try:
            # Process all MSP files
            dfs = []
            for f in input.msps:
                try:
                    df = pd.read_csv(f, sep='\t')
                    dfs.append(df)
                    log_file.write(f"Successfully read {f}, found {len(df)} rows\n")
                except Exception as e:
                    log_file.write(f"Error reading {f}: {str(e)}\n")
            
            if dfs:
                # Combine all dataframes
                combined = pd.concat(dfs, ignore_index=True)
                
                # Sort by chromosome and position
                combined = combined.sort_values(['chm', 'spos', 'epos'])
                
                # Create output directory if it doesn't exist
                os.makedirs(os.path.dirname(output.combined), exist_ok=True)
                
                # Write combined output
                combined.to_csv(output.combined, sep='\t', index=False)
                log_file.write(f"Successfully combined {len(dfs)} MSP files into {output.combined}, total rows: {len(combined)}\n")
            else:
                log_file.write("No valid MSP files found to combine\n")
                # Create empty output file
                with open(output.combined, 'w') as f:
                    f.write("# No valid GnomiX results found to combine\n")
        
        except Exception as e:
            log_file.write(f"Error combining results: {str(e)}\n")
            raise e
        finally:
            log_file.close()

# Generate summary report
rule create_ancestry_report:
    input:
        pca_eigenval="analysis_other/ancestry/pca/merged_cohort.eigenval",
        pca_eigenvec="analysis_other/ancestry/pca/merged_cohort.eigenvec", 
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam",
        iadmix="analysis_other/ancestry/global/iadmix/results.done",
        #admixture="analysis_other/ancestry/global/admixture/results.done",
        rfmix="analysis_other/ancestry/local/rfmix/results.done",
        gnomix="analysis_other/ancestry/local/gnomix/results.done"
    output:
        report="analysis_other/ancestry/reports/ancestry_summary.html"
    conda:
        "../env/r_ancestry.yml"
    log:
        "logs/ancestry/create_ancestry_report.log"
    params:
        rscript = "workflow/scripts/24_generate_ancestry_report.R"
    shell:
        """
        mkdir -p analysis_other/ancestry/reports
        cd analysis_other/ancestry/reports
        Rscript ../../../{params.rscript} >{log} 2>&1
        """
