# Global variables for ancestry analysis
ANCESTRY_TOOLS = ["iadmix", "admixture"]
LOCAL_ANCESTRY_TOOLS = ["rfmix", "gnomix"]
THOUSAND_G_POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

rule all_ancestry:
    input:
        # Merged and QC'd variant set
        "analysis_other/ancestry/plink/merged_cohort.bed",
        "analysis_other/ancestry/plink/merged_cohort.bim", 
        "analysis_other/ancestry/plink/merged_cohort.fam",
        # 1000G reference data
        "analysis_other/ancestry/reference/1000G_phase3_T2T.bed",
        "analysis_other/ancestry/reference/1000G_phase3_T2T.bim",
        "analysis_other/ancestry/reference/1000G_phase3_T2T.fam",
        # Global ancestry results
        expand("analysis_other/ancestry/global/{tool}/results.done", tool=ANCESTRY_TOOLS),
        # Local ancestry results  
        expand("analysis_other/ancestry/local/{tool}/results.done", tool=LOCAL_ANCESTRY_TOOLS),
        # Summary reports
        "analysis_other/ancestry/reports/ancestry_summary.html"

# Download 1000 Genomes Phase 3 data
rule download_1000g_vcf:
    output:
        vcf="data/ref/variant_sets/1000G/chr{chr}_phase3_v5a.vcf.gz",
        tbi="data/ref/variant_sets/1000G/chr{chr}_phase3_v5a.vcf.gz.tbi"
    params:
        url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    shell:
        """
        mkdir -p ancestry/raw_data/1000G
        wget -O {output.vcf} {params.url}
        wget -O {output.tbi} {params.url}.tbi
        """

# Download 1000G sample metadata
rule download_1000g_metadata:
    output:
        metadata="data/ref/variant_sets/1000G/1000G_sample_info.txt",
    shell:
        """
        mkdir -p data/ref/variant_sets/1000G
        
        # Download sample information
        wget -O {output.metadata} \
            http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        """

# Step 1: Add chr prefix and bgzip
rule add_chr_prefix_and_bgzip:
    input:
        vcf="data/ref/variant_sets/1000G/chr{chr}_phase3_v5a.vcf.gz"
    output:
        vcf=temp("ancestry/processed/1000G/chr{chr}_temp.vcf.gz")
    log:
        "logs/ancestry/add_chr_prefix_chr{chr}.log"
    shell:
        """
        zcat {input.vcf} | \
        awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{$1="chr"$1; print}}' | \
        bgzip > {output.vcf} 2>{log}
        """

# Step 2: Liftover coordinates
rule crossmap_liftover:
    input:
        vcf="ancestry/processed/1000G/chr{chr}_temp.vcf.gz"
    output:
        lifted_vcf="data/ref/variant_sets/1000G/chr{chr}_T2T.vcf.gz",
        rejected="data/ref/variant_sets/1000G/chr{chr}_rejected.vcf"
    params:
        chain=config['liftover_chain_GRCh37_to_T2T'],
        ref=config['ref']
    conda:
        "../env/liftover.yml"
    log:
        "logs/ancestry/liftover_1000g_chr{chr}.log"
    shell:
        """
        CrossMap.py vcf {params.chain} \
            {input.vcf} \
            {params.ref} \
            {output.lifted_vcf} >{log} 2>&1
        """

# Step 3: Index lifted VCF
rule index_lifted_vcf:
    input:
        lifted_vcf="data/ref/variant_sets/1000G/chr{chr}_T2T.vcf.gz"
    output:
        tbi="data/ref/variant_sets/1000G/chr{chr}_T2T.vcf.gz.tbi"
    log:
        "logs/ancestry/index_lifted_vcf_chr{chr}.log"
    shell:
        """
        tabix -p vcf {input.lifted_vcf} 2>{log}
        """

# Merge your sample VCFs 
rule merge_sample_vcfs:
    input:
        vcfs=lambda wc: expand("assembly/variants/{sample}/phased_verkko/small_variants.vcf.gz", 
                              sample=[s for s in finished_samples if s in asm_samples])
    output:
        merged="analysis_other/ancestry/processed/samples/merged_samples.vcf.gz",
        list="analysis_other/ancestry/processed/samples/vcf_list.txt"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/ancestry/merge_sample_vcfs.log"
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
            2>{log}
        
        tabix -p vcf {output.merged}
        """

# Convert sample VCFs to Plink format
rule samples_to_plink:
    input:
        vcf="analysis_other/ancestry/processed/samples/merged_samples.vcf.gz"
    output:
        bed="analysis_other/ancestry/plink/samples.bed",
        bim="analysis_other/ancestry/plink/samples.bim", 
        fam="analysis_other/ancestry/plink/samples.fam"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/samples_to_plink.log"
    shell:
        """
        plink2 \
            --vcf {input.vcf} \
            --make-bed \
            --out ancestry/plink/samples \
            --allow-extra-chr \
            2>{log}
        """

# Merge 1000G chromosomes and convert to Plink
rule merge_1000g_to_plink:
    input:
        vcfs=expand("data/ref/variant_sets/1000G/chr{chr}_T2T.vcf.gz", chr=CHROMOSOMES),
        metadata="data/ref/variant_sets/1000G/1000G_sample_info.txt"
    output:
        bed="analysis_other/ancestry/reference/1000G_phase3_T2T.bed",
        bim="analysis_other/ancestry/reference/1000G_phase3_T2T.bim",
        fam="analysis_other/ancestry/reference/1000G_phase3_T2T.fam",
        merged_vcf="analysis_other/ancestry/processed/1000G/1000G_merged_T2T.vcf.gz"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/merge_1000g_to_plink.log"
    shell:
        """
        # Concatenate chromosomes
        bcftools concat \
            {input.vcfs} \
            --output-type z \
            --output {output.merged_vcf} \
            2>{log}
        
        tabix -p vcf {output.merged_vcf}
        
        # Convert to Plink format
        plink2 \
            --vcf {output.merged_vcf} \
            --make-bed \
            --out ancestry/reference/1000G_phase3_T2T \
            --allow-extra-chr \
            2>>{log}
        
        # Update FAM file with population information
        python3 scripts/24_update_1000G.py \
            --fam {output.fam} \
            --metadata {input.metadata} \
            --output {output.fam}.tmp
        
        mv {output.fam}.tmp {output.fam}
        """

# Quality control and harmonization
rule qc_and_harmonize:
    input:
        samples_bed="analysis_other/ancestry/plink/samples.bed",
        samples_bim="analysis_other/ancestry/plink/samples.bim",
        samples_fam="analysis_other/ancestry/plink/samples.fam",
        ref_bed="analysis_other/ancestry/reference/1000G_phase3_T2T.bed",
        ref_bim="analysis_other/ancestry/reference/1000G_phase3_T2T.bim",
        ref_fam="analysis_other/ancestry/reference/1000G_phase3_T2T.fam"
    output:
        bed="analysis_other/ancestry/plink/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged_cohort.fam",
        qc_report="analysis_other/ancestry/qc/harmonization_report.txt"
    conda:
        "../env/plink2.yml"
    log:
        "logs/ancestry/qc_harmonize.log"
    shell:
        """
        # QC filters
        MAF_THRESHOLD=0.01
        GENO_THRESHOLD=0.05
        HWE_THRESHOLD=1e-6
        
        # QC reference data
        plink2 \
            --bfile ancestry/reference/1000G_phase3_T2T \
            --maf $MAF_THRESHOLD \
            --geno $GENO_THRESHOLD \
            --hwe $HWE_THRESHOLD \
            --make-bed \
            --out ancestry/plink/1000G_qc \
            2>{log}
        
        # QC sample data  
        plink2 \
            --bfile ancestry/plink/samples \
            --maf $MAF_THRESHOLD \
            --geno $GENO_THRESHOLD \
            --hwe $HWE_THRESHOLD \
            --make-bed \
            --out ancestry/plink/samples_qc \
            2>>{log}
        
        # Find common variants
        comm -12 \
            <(cut -f2 ancestry/plink/1000G_qc.bim | sort) \
            <(cut -f2 ancestry/plink/samples_qc.bim | sort) \
            > ancestry/plink/common_variants.txt
        
        # Extract common variants from both datasets
        plink2 \
            --bfile ancestry/plink/1000G_qc \
            --extract ancestry/plink/common_variants.txt \
            --make-bed \
            --out ancestry/plink/1000G_common \
            2>>{log}
        
        plink2 \
            --bfile ancestry/plink/samples_qc \
            --extract ancestry/plink/common_variants.txt \
            --make-bed \
            --out ancestry/plink/samples_common \
            2>>{log}
        
        # Merge datasets
        echo "analysis_other/ancestry/plink/samples_common" > ancestry/plink/merge_list.txt
        
        plink2 \
            --bfile ancestry/plink/1000G_common \
            --merge-list ancestry/plink/merge_list.txt \
            --make-bed \
            --out ancestry/plink/merged_cohort \
            2>>{log}
        
        # Generate QC report
        echo "Harmonization completed on $(date)" > {output.qc_report}
        echo "Common variants: $(wc -l < ancestry/plink/common_variants.txt)" >> {output.qc_report}
        echo "Final sample count: $(wc -l < {output.fam})" >> {output.qc_report}
        """

# PCA for population structure
rule run_pca:
    input:
        bed="analysis_other/ancestry/plink/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged_cohort.bim", 
        fam="analysis_other/ancestry/plink/merged_cohort.fam"
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
            --bfile ancestry/plink/merged_cohort \
            --indep-pairwise 50 10 0.2 \
            --out ancestry/pca/ld_pruned \
            2>{log}
        
        # Run PCA
        plink2 \
            --bfile ancestry/plink/merged_cohort \
            --extract ancestry/pca/ld_pruned.prune.in \
            --pca 20 \
            --out ancestry/pca/merged_cohort \
            2>>{log}
        """

# Global ancestry with iAdmix
rule run_iadmix:
    input:
        bed="analysis_other/ancestry/plink/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged_cohort.fam"
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
            --file analysis_other/ancestry/plink/merged_cohort \
            --out analysis_other/ancestry/global/iadmix/results \
            --K 5 \
            2>{log}
        
        # Process results
        cp analysis_other/ancestry/global/iadmix/results.Q {output.results}
        """

# Global ancestry with ADMIXTURE
rule run_admixture:
    input:
        bed="analysis_other/ancestry/plink/merged_cohort.bed"
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
        bed="analysis_other/ancestry/plink/merged_cohort.bed",
        bim="analysis_other/ancestry/plink/merged_cohort.bim",
        fam="analysis_other/ancestry/plink/merged_cohort.fam"
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
            --bfile ancestry/plink/merged_cohort \
            --export vcf bgz \
            --out ancestry/local/input/merged_cohort_phased \
            2>{log}
        
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
            2>{log}
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
            2>{log}
        """

# Generate summary report
rule create_ancestry_report:
    input:
        pca_eigenval="analysis_other/ancestry/pca/merged_cohort.eigenval",
        pca_eigenvec="analysis_other/ancestry/pca/merged_cohort.eigenvec", 
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
    script:
        "../scripts/generate_ancestry_report.R"

