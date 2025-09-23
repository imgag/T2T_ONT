# Global variables for ancestry analysis
ANCESTRY_TOOLS = ["iadmix", "admixture"]
LOCAL_ANCESTRY_TOOLS = ["rfmix", "gnomix"]
THOUSAND_G_POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

rule all_ancestry:
    input:
        # Merged and QC'd variant set
        "analysis_other/ancestry/plink/merged/merged_cohort.bed",
        "analysis_other/ancestry/plink/merged/merged_cohort.bim", 
        "analysis_other/ancestry/plink/merged/merged_cohort.fam",
        # 1000G reference data (filtered)
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.bed",
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.bim",
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.fam",
        "analysis_other/ancestry/plink/reference/1000G_phase3_T2T.afreq",
        # Sample data (filtered)
        "analysis_other/ancestry/plink/samples/samples_filtered.bed",
        "analysis_other/ancestry/plink/samples/samples_filtered.bim",
        "analysis_other/ancestry/plink/samples/samples_filtered.fam",
        # Global ancestry results
        expand("analysis_other/ancestry/global/{tool}/results.done", tool=ANCESTRY_TOOLS),
        # Local ancestry results  
        expand("analysis_other/ancestry/local/{tool}/results.done", tool=LOCAL_ANCESTRY_TOOLS),
        # Summary reports
        "analysis_other/ancestry/reports/ancestry_summary.html"



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

# Step A8: Convert merged VCF to Plink format
rule vcf_to_plink_1000g_with_sex:
    input:
        merged_vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz"
    output:
        bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bed",
        bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bim",
        fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.fam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/ref_vcf_to_plink_1000g_with_sex.log"
    params:
        min_male_xf=config['min_male_xf'],
        max_female_yrate=config['max_female_yrate']
    threads:
        32
    shell:
        """
        plink2 \
            --vcf {input.merged_vcf} \
            --make-bed \
            --out analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --impute-sex max-female-xf={params.max_female_yrate} min-male-xf={params.min_male_xf} \
            --threads {threads} \
            >{log} 2>&1
        """

# Step A9: Caluclate allele frequencies for 1000G data, used in sex imputation of samples
rule generate_1000g_frequencies:
    input:
        bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bed",
        bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bim",
        fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.fam"
    output:
        afreq="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.afreq"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/ref_generate_frequencies.log"
    threads:
        32
    shell:
        """
        plink2 \
            --bfile analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex \
            --freq \
            --out analysis_other/ancestry/plink/reference/1000G_phase3_T2T \
            --threads {threads} \
            >{log} 2>&1
        """



# Step A10: Update FAM file with population information
rule update_1000g_fam:
    input:
        fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.fam",
        metadata="data/ref/variant_sets/1000G/integrated_call_samples_v3.20130502.ALL.panel",
        bed_in="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bed",
        bim_in="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_imputed_sex.bim"
    output:
        fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.fam",
        bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.bed",
        bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.bim"
    log:
        "logs/ancestry/ref_update_fam.log"
    shell:
        """
        python3 workflow/scripts/24_update_1000G.py \
            --fam {input.fam} \
            --metadata {input.metadata} \
            --output {output.fam} \
            2>{log}
        
        # Copy bed and bim files to final location
        cp -v {input.bed_in} {output.bed} >>{log} 2>&1
        cp -v {input.bim_in} {output.bim} >>{log} 2>&1
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

# Step B3: Infer sample with AF from 1000G variants
rule samples_infer_sex_plink:
    input:
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz",
        afreq="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.afreq"

    output:
        bed=temp("analysis_other/ancestry/plink/samples/samples_imputed.bed"),
        bim=temp("analysis_other/ancestry/plink/samples/samples_imputed.bim"),
        fam=temp("analysis_other/ancestry/plink/samples/samples_imputed.fam")
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/samples_infer_sex_plink.log"
    params:
        afreq="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.afreq",
        min_male_xf=config['min_male_xf'],
        max_female_yrate=config['max_female_yrate']
    threads: 32
    shell:
        """
        # Infer sex from X chromosome
        plink2 \
            --vcf {input.vcf} \
            --make-bed \
            --vcf-half-call m \
            --out analysis_other/ancestry/plink/samples/samples_imputed \
            --read-freq {input.afreq} \
            --impute-sex max-female-xf={params.max_female_yrate} min-male-xf={params.min_male_xf} \
            --allow-extra-chr \
            --set-missing-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 60 \
            --split-par b38 \
            --threads {threads} \
            >{log} 2>&1
        """


# Step B4: Update FAM file with pedigree information
rule update_sample_pedigree:
    input:
        bed="analysis_other/ancestry/plink/samples/samples_imputed.bed",
        bim="analysis_other/ancestry/plink/samples/samples_imputed.bim",
        fam="analysis_other/ancestry/plink/samples/samples_imputed.fam"
    output:
        bed="analysis_other/ancestry/plink/samples/samples.bed",
        bim="analysis_other/ancestry/plink/samples/samples.bim", 
        fam="analysis_other/ancestry/plink/samples/samples.fam"
    log:
        "logs/ancestry/sample_update_pedigree.log"
    shell:
        """
        # Update FAM file with pedigree information
        python3 workflow/scripts/27_update_pedigree_info.py \
            --fam {input.fam} \
            --output {output.fam} \
            >{log} 2>&1
        
        # Copy BED and BIM files
        cp {input.bed} {output.bed}
        cp {input.bim} {output.bim}
        """

# C: Merge datasets, QC 

# Step C1: Filter datasets with QC filters. Reassign variant IDs.
rule qc_filter_plink:
    input:
        bed="analysis_other/ancestry/plink/{dataset}/{prefix}.bed",
        bim="analysis_other/ancestry/plink/{dataset}/{prefix}.bim",
        fam="analysis_other/ancestry/plink/{dataset}/{prefix}.fam"
    output:
        bed="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.bed",
        bim="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.bim",
        fam="analysis_other/ancestry/plink/{dataset}/{prefix}_filtered.fam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/qc_filter_{dataset}_{prefix}.log"
    threads: 32
    shell:
        """
        # QC filters - lenient HWE for family structures
        MAF_THRESHOLD=0.01     # Keep rare variants for ancestry
        GENO_THRESHOLD=0.15    # Lenient missing rate
        HWE_THRESHOLD=1e-3     # Lenient HWE for population structure with included families
        
        plink2 \
            --bfile analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix} \
            --maf $MAF_THRESHOLD \
            --geno $GENO_THRESHOLD \
            --hwe $HWE_THRESHOLD \
            --make-bed \
            --out analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix}_filtered \
            --threads {threads} \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 1000 \
            --sort-vars \
            >{log} 2>&1
        """

# Step C2: Merge datasets
rule merge_datasets_with_matching:
    input:
        ref_bed="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.bed",
        ref_bim="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.bim",
        ref_fam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.fam",
        sample_bed="analysis_other/ancestry/plink/samples/samples_filtered.bed",
        sample_bim="analysis_other/ancestry/plink/samples/samples_filtered.bim",
        sample_fam="analysis_other/ancestry/plink/samples/samples_filtered.fam"
    output:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam",
        merge_report="analysis_other/ancestry/qc/merge_report.txt"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/merge_datasets.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/plink/merged analysis_other/ancestry/qc
        
        # Attempt merge with automatic variant matching and strand flipping
        plink2 \
            --bfile analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered \
            --pmerge \
                {input.sample_bed} \
                {input.sample_bim} \
                {input.sample_fam} \
            --make-bed \
            --out analysis_other/ancestry/plink/merged/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        
        # Generate merge report
        echo "Merge completed on $(date)" > {output.merge_report}
        echo "Final sample count: $(wc -l < {output.fam})" >> {output.merge_report}
        echo "Final variant count: $(wc -l < {output.bim})" >> {output.merge_report}
        
        # Check for merge conflicts
        if [[ -f analysis_other/ancestry/plink/merged/merged_cohort.missnp ]]; then
            echo "Variants with merge conflicts: $(wc -l < analysis_other/ancestry/plink/merged/merged_cohort.missnp)" >> {output.merge_report}
        fi
        """

# Step C3: Generate harmonization report
rule generate_harmonization_report:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam",
        common_variants="analysis_other/ancestry/plink/common_variants.txt"
    output:
        qc_report="analysis_other/ancestry/qc/harmonization_report.txt"
    log:
        "logs/ancestry/qc_harmonization_report.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/qc
        
        echo "Harmonization completed on $(date)" > {output.qc_report}
        echo "Common variants: $(wc -l < {input.common_variants})" >> {output.qc_report}
        echo "Final sample count: $(wc -l < {input.fam})" >> {output.qc_report}
        
        # Report family structure
        echo "" >> {output.qc_report}
        echo "Family structure summary:" >> {output.qc_report}
        awk '$3 != "0" || $4 != "0" {{print "Child: " $2 " Father: " $3 " Mother: " $4}}' {input.fam} >> {output.qc_report}
        
        echo "Report generated successfully" >{log}
        """

# Part D: Ancestry and population analyses 

# PCA for population structure (accounting for related samples)
rule run_pca:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim", 
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam"
    output:
        eigenval="analysis_other/ancestry/pca/merged_cohort.eigenval",
        eigenvec="analysis_other/ancestry/pca/merged_cohort.eigenvec"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/run_pca.log"
    shell:
        """
        # LD pruning for PCA
        plink2 \
            --bfile analysis_other/ancestry/plink/merged/merged_cohort \
            --indep-pairwise 50 10 0.2 \
            --out analysis_other/ancestry/pca/ld_pruned \
            >{log} 2>&1
        
        # Run PCA (note: plink2 automatically handles related samples)
        plink2 \
            --bfile analysis_other/ancestry/plink/merged/merged_cohort \
            --extract analysis_other/ancestry/pca/ld_pruned.prune.in \
            --pca 20 \
            --out analysis_other/ancestry/pca/merged_cohort \
            2>>{log}
        """

# Global ancestry with iAdmix
rule run_iadmix:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam"
    output:
        touch("analysis_other/ancestry/global/iadmix/results.done"),
        results="analysis_other/ancestry/global/iadmix/admixture_proportions.txt"
    log:
        "logs/ancestry/iadmix.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/global/iadmix
        
        # Get current working directory and user info
        WORKDIR=$(pwd)
        USER_ID=$(id -u)
        GROUP_ID=$(id -g)
        
        # Run iAdmix using Docker
        docker run --rm \
            --user ${{USER_ID}}:${{GROUP_ID}} \
            -v ${{WORKDIR}}:/workdir \
            -w /workdir \
            caspargross/iadmix \
            iadmix \
            --file analysis_other/ancestry/plink/merged/merged_cohort \
            --out analysis_other/ancestry/global/iadmix/results \
            --K 5 \
            >{log} 2>&1
        
        # Process results
        cp analysis_other/ancestry/global/iadmix/results.Q {output.results}
        """

# Global ancestry with ADMIXTURE
rule run_admixture:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed"
    output:
        touch("analysis_other/ancestry/global/admixture/results.done"),
        results="analysis_other/ancestry/global/admixture/merged_cohort.5.Q"
    conda:
        "../env/admixture.yml"
    log:
        "logs/ancestry/admixture.log"
    threads: 4
    shell:
        """
        cd ancestry/global/admixture
        
        # Run ADMIXTURE for different K values
        for K in {{3..7}}; do
            admixture --cv ../../../{input.bed} $K -j{threads} 2>>../../../{log}
        done
        
        # The output files will be created in the working directory
        """

# Prepare data for local ancestry analysis
rule prepare_local_ancestry:
    input:
        bed="analysis_other/ancestry/plink/merged/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam"
    output:
        phased_vcf="analysis_other/ancestry/local/input/merged_cohort_phased.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/prepare_local_ancestry.log"
    shell:
        """
        # Convert back to VCF format (assuming phased data)
        plink2 \
            --bfile ancestry/plink/merged/merged_cohort \
            --export vcf bgz \
            --out ancestry/local/input/merged_cohort_phased \
            >{log} 2>&1
        
        # Create sample mapping file for RFMix
        python3 scripts/25_create_ancestry_sample_map.py \
            --fam {input.fam} \
            --output {output.sample_map}
        """

# Local ancestry with RFMix
rule run_rfmix:
    input:
        vcf="analysis_other/ancestry/local/input/merged_cohort_phased.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt"
    output:
        touch("analysis_other/ancestry/local/rfmix/results.done"),
        results="analysis_other/ancestry/local/rfmix/output.rfmix.Q"
    conda:
        "../env/rfmix.yml"
    log:
        "logs/ancestry/rfmix.log"
    threads: 8
    shell:
        """
        mkdir -p ancestry/local/rfmix
        
        # Run RFMix
        rfmix \
            -f {input.vcf} \
            -r {input.sample_map} \
            -m 5 \
            -g ancestry/local/rfmix/genetic_map.txt \
            -o ancestry/local/rfmix/output \
            --chromosome=ALL \
            --n-threads={threads} \
            >{log} 2>&1
        """

# Local ancestry with gnomix  
rule run_gnomix:
    input:
        vcf="analysis_other/ancestry/local/input/merged_cohort_phased.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt"
    output:
        touch("analysis_other/ancestry/local/gnomix/results.done"),
        results="analysis_other/ancestry/local/gnomix/output.lai"
    conda:
        "../env/gnomix.yml"
    log:
        "logs/ancestry/gnomix.log"
    threads: 8
    shell:
        """
        mkdir -p ancestry/local/gnomix
        
        # Run gnomix
        python gnomix.py \
            --query_file {input.vcf} \
            --output_basename ancestry/local/gnomix/output \
            --population_map {input.sample_map} \
            --n_cores {threads} \
            >{log} 2>&1
        """

# Generate summary report
rule create_ancestry_report:
    input:
        pca_eigenval="analysis_other/ancestry/pca/merged_cohort.eigenval",
        pca_eigenvec="analysis_other/ancestry/pca/merged_cohort.eigenvec", 
        fam="analysis_other/ancestry/plink/merged/merged_cohort.fam",
        iadmix="analysis_other/ancestry/global/iadmix/results.done",
        admixture="analysis_other/ancestry/global/admixture/results.done",
        rfmix="analysis_other/ancestry/local/rfmix/results.done",
        gnomix="analysis_other/ancestry/local/gnomix/results.done"
    output:
        report="analysis_other/ancestry/reports/ancestry_summary.html"
    conda:
        "../env/r_ancestry.yml"
    log:
        "logs/ancestry/create_ancestry_report.log"
    params:
        rscript = "../../../scripts/26_generate_ancestry_report.R"
    shell:
        """
        cd analysis_other/ancestry/reports
        Rscript {params.rscript} >{log} 2>&1
        """
