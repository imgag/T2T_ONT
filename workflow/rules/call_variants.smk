# Determine sample sex using XY coverage ratio
rule determine_sex:
    input:
        bam=lambda wc: re.sub(
            ".chr\d+",
            "",
            get_assembly_input(wc)["hq"]
            .replace("assembly/input/", "data/mapped/")
            .replace(".fastq.gz", ".bam"),
        ),
    output:
        sex="assembly/qc/{isphased}_{tool}/{asm}/sample_sex.txt",
    log:
        "logs/determine_sex_{isphased}_{tool}_{asm}.log",
    conda:
        "../env/mosdepth.yml"
    shell:
        """
        # Calculate X chromosome coverage
        mosdepth -n --by 10000 --fast-mode --chrom chrX \
            $(dirname {output.sex})/$(basename {output.sex}.tmp.X) {input.bam} 2>{log}

        # Calculate Y chromosome coverage 
        mosdepth -n --by 10000 --fast-mode --chrom chrY \
            $(dirname {output.sex})/$(basename {output.sex}.tmp.Y) {input.bam} 2>>{log}

        # Get mean coverage for X and Y
        X_cov=$(zcat $(dirname {output.sex})/$(basename {output.sex}.tmp.X.regions.bed.gz) | \
            awk '{{sum+=$4}} END {{print sum/NR}}')
        Y_cov=$(zcat $(dirname {output.sex})/$(basename {output.sex}.tmp.Y.regions.bed.gz) | \
            awk '{{sum+=$4}} END {{print sum/NR}}')

        # Determine sex based on Y/X ratio
        # Ratio > 0.5 indicates male, < 0.2 indicates female
        ratio=$(echo "scale=3; $Y_cov/$X_cov" | bc)
        if (( $(echo "$ratio > 0.5" | bc -l) )); then
            echo "male\t${{ratio}}" > {output.sex}
        else
            echo "female\t${{ratio}}" > {output.sex}
        fi

        # Cleanup temp files
        rm $(dirname {output.sex})/$(basename {output.sex}.tmp)*
        """


def get_sex_file(wildcards):
    # Always use "phased" since dipcall is only used for phased assemblies
    sex_file = f"assembly/qc/phased_{wildcards.tool}/{wildcards.asm}/sample_sex.txt"
    return (
        "-x " + config["X_PAR_file"] if open(sex_file).read().strip() == "male" else ""
    )


# Small Variants and Indels with dipcall
rule dipcall:
    input:
        pat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype1"}),
        mat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype2"}),
        ref_fa=get_ref_genome,
        sex="assembly/qc/phased_{tool}/{asm}/sample_sex.txt",
    output:
        makefile="assembly/variants/{asm}/phased_{tool}/small_variants.dip.mak",
        vcf="assembly/variants/{asm}/phased_{tool}/small_variants.dip.vcf.gz",
    params:
        sex=get_sex_file,
        run_dipcall=config["run-dipcall"],
    log:
        "logs/dipcall_{tool}_{asm}.log",
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


# Structural Variants with hapdiff (minimap2+svim_asm)
rule hapdiff:
    input:
        ref_fa=get_ref_genome,
        pat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype1"}),
        mat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype2"}),
    output:
        vcf="assembly/variants/{asm}/phased_{tool}/hapdiff_phased.vcf.gz",
    log:
        "logs/hapdiff_{tool}_{asm}.log",
    conda:
        "../env/hapdiff.yml"
    threads: 20
    params:
        hapdiff=config["hapdiff"],
    shell:
        """
        {params.hapdiff} \
            --reference {input.ref_fa} \
            --pat {input.pat_fa} \
            --mat {input.mat_fa} \
            --out-dir $(dirname {output.vcf}) \
            -t {threads} \
            >{log} 2>&1
        """


# whatshap stats
rule whatshap_stats:
    input:
        vcf=rules.dipcall.output.vcf,
    output:
        stats="assembly/qc/phased_{tool}/{asm}/whatshap_stats.txt",
        tsv="assembly/qc/phased_{tool}/{asm}/whatshap_stats.tsv",
        blocks="assembly/qc/phased_{tool}/{asm}/whatshap_blocks.tsv",
        gtf="assembly/qc/phased_{tool}/{asm}/whatshap_blocks.gtf",
    log:
        "logs/whatshap_stats_{tool}_{asm}.log",
    conda:
        "../env/whatshap.yml"
    shell:
        """
        whatshap stats \
            --tsv={output.tsv} \
            --block-list={output.blocks} \
            --gtf={output.gtf} \
            {input.vcf} \
            >{output.stats} 2>&1
        """


# Filter VCFs to only include variants present in both samples
rule filter_shared_variants:
    input:
        ref=config["ref_giab_vcf"],
        vcf=rules.dipcall.output.vcf,
    output:
        ref_filtered="assembly/qc/phased_{tool}/{asm}/giab.filtered.vcf.gz",
        asm_filtered="assembly/qc/phased_{tool}/{asm}/assembly.filtered.vcf.gz",
    log:
        "logs/filter_shared_variants_{tool}_{asm}.log",
    conda:
        "../env/bcftools.yml"
    shell:
        """
        # Extract shared variants and normalize sample names
        temp_ref=$(mktemp)
        temp_asm=$(mktemp)

        bcftools view -t ^chrX,^chrY {input.ref} | \
        bcftools reheader --samples <(echo "SAMPLE") | \
        bgzip > "$temp_ref" && \
        bcftools index "$temp_ref"

        bcftools view -t ^chrX,^chrY {input.vcf} | \
        bcftools reheader --samples <(echo "SAMPLE") | \
        bgzip > "$temp_asm" && \
        bcftools index "$temp_asm"

        # Get intersection of variants
        bcftools isec -n=2 -w1 "$temp_ref" "$temp_asm" | \
        bcftools view -Oz > {output.ref_filtered}

        bcftools isec -n=2 -w1 "$temp_asm" "$temp_ref" | \
        bcftools view -Oz > {output.asm_filtered}

        # Clean up temp files
        rm "$temp_ref" "$temp_asm"

        # Index output files
        bcftools index {output.ref_filtered}
        bcftools index {output.asm_filtered}
        """


# whatshap compare for phasing accuracy
rule whatshap_compare:
    input:
        ref=rules.filter_shared_variants.output.ref_filtered,
        vcf=rules.filter_shared_variants.output.asm_filtered,
    output:
        stats="assembly/qc/phased_{tool}/{asm}/whatshap_compare.txt",
        tsv="assembly/qc/phased_{tool}/{asm}/whatshap_compare.tsv",
    log:
        "logs/whatshap_compare_{tool}_{asm}.log",
    conda:
        "../env/whatshap.yml"
    shell:
        """
        whatshap compare \
            --names giab,assembly \
            --tsv-pairwise {output.tsv} \
            {input.ref} \
            {input.vcf} \
            >{output.stats} 2>{log}
        """
