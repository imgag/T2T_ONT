# Global variables for ancestry analysis
ANCESTRY_TOOLS = ["iadmix", "admixture"]
LOCAL_ANCESTRY_TOOLS = ["rfmix", "gnomix"]
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

# Step A8: Convert VCF to PGEN format
rule vcf_to_plink_1000g:
    input:
        merged_vcf="analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz"
    output:
        pgen=temp("analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.pgen"),
        pvar=temp("analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.pvar"),
        psam=temp("analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.psam")
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
            --out analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --threads {threads} \
            >{log} 2>&1
        """

# Step A9: Update 1000G PSAM with metadata  
rule update_1000g_psam:
    input:
        pgen="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.pgen",
        pvar="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.pvar",
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_temp.psam",
        metadata="data/ref/variant_sets/1000G/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        pgen="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.pgen",
        pvar="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.pvar",
        psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T.psam"
    log:
        "logs/ancestry/update_1000g_psam.log"
    shell:
        """
        # Convert 1000G metadata to PSAM format
        awk 'BEGIN{{OFS="\\t"; print "#FID","IID","PAT","MAT","SEX","PHENO"}} 
             NR>1 {{print $1,$2,$3,$4,$5,$6}}' \
             {input.metadata} > {output.psam} 2>{log}
        
        # Copy pgen and pvar files
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
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
        vcf="analysis_other/ancestry/processed_vcf/samples/merged_samples_filtered.vcf.gz"
    output:
        pgen=temp("analysis_other/ancestry/plink/samples/samples_temp.pgen"),
        pvar=temp("analysis_other/ancestry/plink/samples/samples_temp.pvar"),
        psam=temp("analysis_other/ancestry/plink/samples/samples_temp.psam")
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
            --out analysis_other/ancestry/plink/samples/samples_temp \
            --allow-extra-chr \
            --split-par b38 \
            --vcf-half-call m \
            --set-missing-var-ids @:# \
            --threads {threads} \
            >{log} 2>&1
        """

# Step B4: Update sample PSAM with metadata
rule update_sample_psam:
    input:
        pgen="analysis_other/ancestry/plink/samples/samples_temp.pgen",
        pvar="analysis_other/ancestry/plink/samples/samples_temp.pvar", 
        psam="analysis_other/ancestry/plink/samples/samples_temp.psam",
        metadata="data/samples_pedigree.tsv"
    output:
        pgen="analysis_other/ancestry/plink/samples/samples.pgen",
        pvar="analysis_other/ancestry/plink/samples/samples.pvar",
        psam="analysis_other/ancestry/plink/samples/samples.psam"
    log:
        "logs/ancestry/update_sample_psam.log"
    shell:
        """
        # Use provided metadata for sample pedigree
        cp {input.metadata} {output.psam} 2>{log}
        
        # Copy pgen and pvar files  
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
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
    shell:
        """
        # QC filters - lenient for ancestry analysis
        MAF_THRESHOLD=0.01     # Keep rare variants for ancestry
        GENO_THRESHOLD=0.15    # Lenient missing rate
        HWE_THRESHOLD=1e-3     # Lenient HWE for population structure with families
        MIND_THRESHOLD=0.1     # Remove samples with >10% missing data
        
        echo "Applying QC filters to {wildcards.dataset}/{wildcards.prefix}..." >{log}
        echo "GENO_THRESHOLD=$GENO_THRESHOLD (allowing up to 15% missing per variant)" >>{log}
        
        plink2 \
            --pfile analysis_other/ancestry/plink/{wildcards.dataset}/{wildcards.prefix} \
            --mind $MIND_THRESHOLD \
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

# Step C2: Merge datasets using PGEN format
rule merge_datasets_with_matching:
    input:
        ref_pgen="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.pgen",
        ref_pvar="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.pvar",
        ref_psam="analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered.psam",
        sample_pgen="analysis_other/ancestry/plink/samples/samples_filtered.pgen",
        sample_pvar="analysis_other/ancestry/plink/samples/samples_filtered.pvar",
        sample_psam="analysis_other/ancestry/plink/samples/samples_filtered.psam"
    output:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam",
        merge_report="analysis_other/ancestry/qc/merge_report.txt"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/merge_datasets.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/plink/merged analysis_other/ancestry/qc
        
        # Create pmerge list file
        echo "analysis_other/ancestry/plink/samples/samples_filtered" > analysis_other/ancestry/plink/pmerge_list.txt
        
        # Attempt merge with automatic variant matching
        plink2 \
            --pfile analysis_other/ancestry/plink/reference/1000G_phase3_T2T_filtered \
            --pmerge-list analysis_other/ancestry/plink/pmerge_list.txt \
            --pmerge-list-mode 6 \
            --make-pgen \
            --out analysis_other/ancestry/plink/merged/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        
        # Generate merge report
        echo "Merge completed on $(date)" > {output.merge_report}
        echo "Final sample count: $(tail -n +2 {output.psam} | wc -l)" >> {output.merge_report}
        echo "Final variant count: $(tail -n +2 {output.pvar} | wc -l)" >> {output.merge_report}
        
        # Check for merge conflicts
        if [[ -f analysis_other/ancestry/plink/merged/merged_cohort.missnp ]]; then
            echo "Variants with merge conflicts: $(wc -l < analysis_other/ancestry/plink/merged/merged_cohort.missnp)" >> {output.merge_report}
        fi
        """

# Step C3: Generate harmonization report (updated for PGEN format)
rule generate_harmonization_report:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        qc_report="analysis_other/ancestry/qc/harmonization_report.txt"
    log:
        "logs/ancestry/qc_harmonization_report.log"
    shell:
        """
        mkdir -p analysis_other/ancestry/qc
        
        echo "Harmonization completed on $(date)" > {output.qc_report}
        echo "Final sample count: $(tail -n +2 {input.psam} | wc -l)" >> {output.qc_report}
        echo "Final variant count: $(tail -n +2 {input.pvar} | wc -l)" >> {output.qc_report}
        
        # Report family structure from PSAM file
        echo "" >> {output.qc_report}
        echo "Family structure summary:" >> {output.qc_report}
        awk 'NR>1 && ($3 != "0" || $4 != "0") {{print "Child: " $2 " Father: " $3 " Mother: " $4}}' {input.psam} >> {output.qc_report}
        
        echo "Report generated successfully" >{log}
        """

# PCA for population structure using PGEN format
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
            2>>{log}
        """

# Global ancestry with iAdmix (convert to BED format for compatibility)
rule run_iadmix:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        touch("analysis_other/ancestry/global/iadmix/results.done"),
        results="analysis_other/ancestry/global/iadmix/admixture_proportions.txt",
        bed=temp("analysis_other/ancestry/global/iadmix/merged_cohort.bed"),
        bim=temp("analysis_other/ancestry/global/iadmix/merged_cohort.bim"),
        fam=temp("analysis_other/ancestry/global/iadmix/merged_cohort.fam")
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/iadmix.log"
    threads: 32
    shell:
        """
        mkdir -p analysis_other/ancestry/global/iadmix
        
        # Convert PGEN to BED format for iAdmix compatibility
        plink2 \
            --pfile analysis_other/ancestry/plink/merged/merged_cohort \
            --make-bed \
            --out analysis_other/ancestry/global/iadmix/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        
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
            --file analysis_other/ancestry/global/iadmix/merged_cohort \
            --out analysis_other/ancestry/global/iadmix/results \
            --K 5 \
            2>>{log}
        
        # Process results
        cp analysis_other/ancestry/global/iadmix/results.Q {output.results} 2>>{log} || echo "Results.Q not found" >>{log}
        """

# Global ancestry with ADMIXTURE (convert to BED format for compatibility)
rule run_admixture:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        touch("analysis_other/ancestry/global/admixture/results.done"),
        results="analysis_other/ancestry/global/admixture/merged_cohort.5.Q",
        bed=temp("analysis_other/ancestry/global/admixture/merged_cohort.bed")
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/admixture.log"
    threads: 4
    shell:
        """
        mkdir -p analysis_other/ancestry/global/admixture
        
        # Convert PGEN to BED format for ADMIXTURE compatibility
        plink2 \
            --pfile analysis_other/ancestry/plink/merged/merged_cohort \
            --make-bed \
            --out analysis_other/ancestry/global/admixture/merged_cohort \
            --threads {threads} \
            >{log} 2>&1
        
        cd analysis_other/ancestry/global/admixture
        
        # Run ADMIXTURE for different K values
        for K in {{3..7}}; do
            admixture --cv merged_cohort.bed $K -j{threads} 2>>../../../../{log}
        done
        """

# Prepare data for local ancestry analysis (convert from PGEN to VCF)
rule prepare_local_ancestry:
    input:
        pgen="analysis_other/ancestry/plink/merged/merged_cohort.pgen",
        pvar="analysis_other/ancestry/plink/merged/merged_cohort.pvar",
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam"
    output:
        phased_vcf="analysis_other/ancestry/local/input/merged_cohort_phased.vcf.gz",
        sample_map="analysis_other/ancestry/local/input/sample_map.txt"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/prepare_local_ancestry.log"
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
        
        # Create sample mapping file for RFMix
        python3 workflow/scripts/25_create_ancestry_sample_map.py \
            --psam {input.psam} \
            --output {output.sample_map} \
            2>>{log}
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
        psam="analysis_other/ancestry/plink/merged/merged_cohort.psam",
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
        rscript = "workflow/scripts/24_generate_ancestry_report.R"
    shell:
        """
        mkdir -p analysis_other/ancestry/reports
        cd analysis_other/ancestry/reports
        Rscript ../../../{params.rscript} >{log} 2>&1
        """
