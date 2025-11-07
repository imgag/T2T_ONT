def get_chr_list_for_asm(wc):

    file = f"assembly/qc/phased_verkko/{wc.asm}/sample_sex.txt"

    if open(file).read().strip() == "male":
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X", "Y"]]
    else:
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X"]]

    return(chroms)

rule all_dip3d:
    input:
        expand("analysis_other/dip3d/{asm}/dip3d.done", asm = finished_samples),
        "analysis_other/dip3d/combined_dip3d_stats.tsv"

rule combine_dip3d_stats:
    input:
        imputed_stats=expand("analysis_other/dip3d/{asm}/4-haplotag/imputed_dip3d_stats.txt", asm=finished_samples),
        snp_tagged_stats=expand("analysis_other/dip3d/{asm}/4-haplotag/snp-tagged_dip3d_stats.txt", asm=finished_samples)
    output:
        combined="analysis_other/dip3d/combined_dip3d_stats.tsv"
    conda:
        "../env/py_report.yml"
    log:
        "logs/dip3d/combine_stats/combine_dip3d_stats.log"
    shell:
        """
        python3 workflow/scripts/42_parse_dip3d_stats.py \
            --imputed {input.imputed_stats} \
            --snp-tagged {input.snp_tagged_stats} \
            --output {output.combined} \
            >{log} 2>&1
        """

rule collect_dip3d:
    input:
        "assembly/qc/phased_verkko/{asm}/sample_sex.txt",
        expand("analysis_other/dip3d/{{asm}}/4-haplotag/{type}_dip3d_stats.txt", type = ["imputed", "snp-tagged"]),
        expand("analysis_other/dip3d/{{asm}}/4-haplotag/{{asm}}.{hp}.bam", hp = ["hp1", "hp2", "tagged", "tagged.sorted"])
    #   lambda wc: expand("analysis_other/dip3d/{{asm}}/5-ashic/{chr}/ashic_result.txt", chr = get_chr_list_for_asm(wc))
    output:
        "analysis_other/dip3d/{asm}/dip3d.done"
    benchmark:
        "runtimes/dip3d/all_merge_dip3d/{asm}.all_merge_dip3d.txt"
    shell:
        """
        touch {output}
        """

# Step 1: Mapping
# ----------------
# Map using Falign (Dip3d paper)
# Then split into chromosomes, and create subset passing dip3d filter criteria

rule dip3d_build_repeat_regions:
    input:
        ref = config['ref']
    output:
        bed = config['ref'] + ".repeat_regions.bed"
    params:
        falign = config['falign']
    log:
        "logs/dip3d/build_repeat_regions/build_repeat_regions.log"
    benchmark:
        "runtimes/dip3d/dip3d_build_repeat_regions/build_repeat_regions.dip3d_build_repeat_regions.txt"
    shell:
        """
        {params.falign} build-repeat \
            {input.ref} \
            {output.bed} \
            >{log} 2>&1
        """

rule dip3d_map_reads:
    input:
        fq = lambda wc: f"assembly/input/{asm[wc.asm]["dataset"]}/{asm[wc.asm]["dataset"]}.POREC_all.fastq.gz",
        fa = config['ref'],
        bed = config['ref'] + ".repeat_regions.bed"
    output:
        bam = temp("analysis_other/dip3d/{asm}/1-falign/porec.fragments.bam")
    params:
        falign = config['falign']
    threads: 30
    log:
        "logs/dip3d/map_reads/{asm}.log"
    benchmark:
        "runtimes/dip3d/dip3d_map_reads/{asm}.dip3d_map_reads.txt"
    shell:
        """
        {params.falign} \
            -repeat_bed {input.bed} \
            -num_threads {threads} \
            -outfmt frag-bam \
            -out {output.bam} \
            {input.fa} \
            {input.fq} \
            >{log} 2>&1
        """

# Step 2:
# -------------------------------
# Split Bam into chromosomes, calculate coverage, apply quality filters


rule dip3d_split_bams:
    input:
        bam = rules.dip3d_map_reads.output.bam
    output:
        temp("analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam")
    params:
        dip3d = config['dip3d'],
    threads: 16
    log:
        "logs/dip3d/split_bams/{asm}.{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_split_bams/{asm}.{chr}.dip3d_split_bams.txt"
    shell:
        """
        {params.dip3d} split-bam \
            {input.bam} \
            analysis_other/dip3d/{wildcards.asm}/2-chr-bams \
            {wildcards.chr} \
            >{log} 2>&1
        """

rule dip3d_select_chr_snp_bam:
    input:
        bam_dirs = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam"
    output:
        bam = temp("analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.hq.bam")
    params:
        dip3d = config['dip3d'],
        l = 100,
        q = 5,
        i = 90,
        cov = 5
    threads: 16
    log:
        "logs/dip3d/select_chr_snp_bam/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_select_chr_snp_bam/{asm}.{chr}.dip3d_select_chr_snp_bam.txt"
    shell:
        """
        {params.dip3d} select-chr-snp-bam \
            -l {params.l} \
            -q {params.q} \
            -i {params.i} \
            analysis_other/dip3d/{wildcards.asm}/2-chr-bams \
            {params.cov} \
            {output.bam} \
            {wildcards.chr} \
            >{log} 2>&1
        """

rule dip3d_calculate_coverage:
    input:
        bam = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.hq.bam"
    output:
        mosdepth_bed = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.per-base.bed.gz"
    params:
        prefix = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}",
        min_coverage = 5
    conda:
        "../env/mosdepth.yml"
    threads: 4
    log:
        "logs/dip3d/calculate_coverage/{asm}.{chr}.calculate_coverage.log"
    benchmark:
        "runtimes/dip3d/dip3d_calculate_coverage/{asm}.{chr}.dip3d_calculate_coverage.txt"
    shell:
        """
        # Calculate per-base coverage with mosdepth
        mosdepth \
            -t {threads} \
            {params.prefix} {input.bam} \
            >{log} 2>&1
        """

rule dip3d_create_well_covered_regions:
    input:
        coverage_bed = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.per-base.bed.gz"
    output:
        bed = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.well_covered.bed"
    params:
        min_coverage = 5
    log:
        "logs/dip3d/create_well_covered_regions/{asm}.{chr}.create_well_covered_regions.log"
    benchmark:
        "runtimes/dip3d/dip3d_create_well_covered_regions/{asm}.{chr}.dip3d_create_well_covered_regions.txt"
    run:
        import gzip

        def create_well_covered_regions(coverage_file, output_file, min_coverage, log_file):
            with open(log_file, 'w') as log:
                log.write(f"Processing coverage file: {coverage_file}\n")
                log.write(f"Minimum coverage threshold: {min_coverage}\n")

                regions = []
                current_region = None

                with gzip.open(coverage_file, 'rt') as f:
                    for line_num, line in enumerate(f, 1):
                        if line_num % 1000000 == 0:
                            log.write(f"Processed {line_num} lines\n")
                            log.flush()

                        parts = line.strip().split('\t')
                        if len(parts) < 4:
                            continue

                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        coverage = float(parts[3])

                        if coverage >= min_coverage:
                            if current_region is None:
                                # Start new region
                                current_region = [chrom, start, end]
                            elif (current_region[0] == chrom and
                                  current_region[2] == start):
                                # Extend current region
                                current_region[2] = end
                            else:
                                # Save current region and start new one
                                regions.append(current_region)
                                current_region = [chrom, start, end]
                        else:
                            # Coverage below threshold, close current region if exists
                            if current_region is not None:
                                regions.append(current_region)
                                current_region = None

                # Don't forget the last region
                if current_region is not None:
                    regions.append(current_region)

                log.write(f"Found {len(regions)} well-covered regions\n")

                # Write output BED file
                with open(output_file, 'w') as out:
                    for region in regions:
                        out.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

                log.write(f"Output written to: {output_file}\n")

        # Execute the function
        create_well_covered_regions(
            input.coverage_bed,
            output.bed,
            params.min_coverage,
            log[0]
        )


# Step 3: Process SNP Set
# -------------------------------------------------------------------
# Split SNP set to chrom, filter to contain only regions well covered by reads
# Convert phasing format to PS tag using whatshap

rule dip3d_split_vcf:
    input:
        vcf = "assembly/variants/{asm}/phased_verkko/small_variants.dip.vcf.gz"
    output:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.vcf"
    params:
        dip3d = config['dip3d'],
    log:
        "logs/dip3d/split_vcf/{asm}.{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_split_vcf/{asm}.{chr}.dip3d_split_vcf.txt"
    shell:
        """
        {params.dip3d} split-vcf \
            -v -p \
            {input.vcf} \
            analysis_other/dip3d/{wildcards.asm}/3-snp \
            {wildcards.chr} \
            >{log} 2>&1
        """

rule dip3d_filter_vcf_by_coverage:
    input:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.vcf",
        bed = rules.dip3d_create_well_covered_regions.output.bed
    output:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.filtered.vcf"
    conda:
        "../env/bcftools.yml"
    log:
        "logs/dip3d/filter_vcf_by_coverage/{asm}.{chr}.filter_vcf_by_coverage.log"
    benchmark:
        "runtimes/dip3d/dip3d_filter_vcf_by_coverage/{asm}.{chr}.dip3d_filter_vcf_by_coverage.txt"
    shell:
        """
        # Compress and index VCF file first
        bgzip -c {input.vcf} > {input.vcf}.gz 2>>{log}
        tabix -p vcf {input.vcf}.gz 2>>{log}

        # Intersect VCF with coverage BED file to keep only well-covered regions
        bcftools view -R {input.bed} {input.vcf}.gz > {output.vcf} 2>>{log}

        # Clean up temporary files
        rm -f {input.vcf}.gz {input.vcf}.gz.tbi 2>>{log}
        """

rule dip3d_whatshap_reformat:
    input:
        ref = config['ref'],
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.filtered.vcf"
    output:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.whatshap.vcf"
    conda:
        "../env/whatshap.yml"
    log:
        "logs/dip3d/whatshap_phase_vcf/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/whatshap_phase_vcf/{asm}.{chr}.whatshap_phase_vcf.txt"
    shell:
        """
        whatshap phase \
            -o {output.vcf} \
            --reference {input.ref} \
            {input.vcf} {input.vcf} \
            >{log} 2>&1
        """

# Step 4
# -------------------
# Haplotag bam file

rule dip3d_bam_haplotag:
    input:
        ref = config['ref'],
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.whatshap.vcf",
        bam = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam"
    output:
        bam = temp("analysis_other/dip3d/{asm}/4-haplotag/{chr}/tagged.bam"),
        imputed_frag_list = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/imputed-frag-hap-list",
        tagged_frag_list = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/snp-tagged-frag-hap-list"
    params:
        dip3d = config['dip3d']
    threads: 8
    log:
        "logs/dip3d/bam_haplotag/{asm}.{chr}.bam_haplotag.log"
    benchmark:
        "runtimes/dip3d/dip3d_bam_haplotag/{asm}.{chr}.haplotag.txt"
    shell:
        """
        {params.dip3d} tag-bam \
            $(dirname {output.bam}) \
            {input.ref} \
            {input.vcf} \
            {input.bam} \
            -t {threads} \
            >{log} 2>&1
        """

rule dip3d_stats_all:
    input:
        lambda wc: expand("analysis_other/dip3d/{{asm}}/4-haplotag/{chr}/{{type}}-frag-hap-list", chr = get_chr_list_for_asm(wc))
    output:
        "analysis_other/dip3d/{asm}/4-haplotag/{type}_dip3d_stats.txt"
    params:
        dip3d = config['dip3d']
    threads: 1
    log:
        "logs/dip3d/stats/{asm}.{type}.stats.log"
    benchmark:
        "runtimes/dip3d/dip3d_bam_haplotag/{asm}.{type}.haplotag.txt"
    shell:
        """
        {params.dip3d} haplo-tag-stats \
            {input} \
            2>{output}
        """

rule dip3d_extract_hp_bams:
    input:
        bam = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/tagged.bam"
    output:
        bam = temp("analysis_other/dip3d/{asm}/4-haplotag/{chr}/hp{hp}.bam"),
    conda:
        "../env/samtools.yml"
    log:
        "logs/dip3d/4-haplotag/{asm}.{chr}.extract_hp{hp}_bams.log"
    shell:
        '''
        samtools view -h {input.bam} | awk '($0 ~ /^@/ || $0 ~ /HP:i:{wildcards.hp}/) {{print $0}}' | samtools view -b -o {output.bam} - 2>>{log}
        '''

rule dip3d_merge_bam:
    input:
        lambda wc: expand(
            "analysis_other/dip3d/{asm}/4-haplotag/{chr}/tagged.bam" if wc.hp == "tagged" else "analysis_other/dip3d/{asm}/4-haplotag/{chr}/{hp}.bam" ,
            asm=wc.asm,
            hp=wc.hp,
            chr=get_chr_list_for_asm(wc)
        )
    output:
        "analysis_other/dip3d/{asm}/4-haplotag/{asm}.{hp}.bam"
    threads: 8
    log:
        "logs/dip3d/merge_bam/{asm}.{hp}.merge_bam.log"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools merge -O BAM - {input} > {output}
        """

rule dip3d_sort_bam:
    input:
        "analysis_other/dip3d/{asm}/4-haplotag/{asm}.{hp}.bam"
    output:
        "analysis_other/dip3d/{asm}/4-haplotag/{asm}.{hp}.sorted.bam"
    threads: 8
    log:
        "logs/dip3d/merge_bam/{asm}.{hp}.sort_bam.log"
    conda:
        "../env/samtools.yml"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2>{log}
        samtools index {output} 2>>{log}
        """


## Step 5: Chromatin contact modeling (ASHIC)
# --------------------------------------

rule dip3d_get_chr_sizes:
    input:
        ref = config['ref']
    output:
        txt = config['ref'] + ".chrom-size.txt"
    params:
        dip3d = config['dip3d']
    log:
        "logs/dip3d/get_chr_sizes/get_chr_sizes.log"
    benchmark:
        "runtimes/dip3d/dip3d_get_chr_sizes/get_chr_sizes.dip3d_get_chr_sizes.txt"
    shell:
        """
        {params.dip3d} extract-ashic-chr-size \
            {input.ref} \
            {output.txt} \
            >{log} 2>&1
        """

rule dip3d_frag_to_ashic_read_pair:
    input:
        frag_list = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/imputed-frag-hap-list"
    output:
        pairs = "analysis_other/dip3d/{asm}/5-ashic/{chr}/ashic_read_pair"
    params:
        dip3d = config['dip3d'],
        min_mapQ = 5,
        chr = lambda wc: wc.chr
    log:
        "logs/dip3d/frag_to_ashic_read_pair/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_frag_to_ashic_read_pair/{asm}.{chr}.dip3d_frag_to_ashic_read_pair.txt"
    shell:
        """
        {params.dip3d} frag-to-ashic-read-pair \
            {input.frag_list} \
            {params.chr} \
            {params.min_mapQ} \
            {output} \
            >{log} 2>&1
        """

rule dip3d_ashic_split2chrs:
    input:
        read_pair = rules.dip3d_frag_to_ashic_read_pair.output
    output:
        directory("analysis_other/dip3d/{asm}/5-ashic/{chr}/split")
    conda:
        "../env/ashic.yml"
    params:
        ashic = config['ashic_data'],
        chr = lambda wc: wc.chr
    log:
        "logs/dip3d/ashic_split2chrs/{asm}.{chr}.split2chrs.log"
    benchmark:
        "runtimes/dip3d/dip3d_ashic_split2chrs/{asm}.{chr}.dip3d_ashic_split2chrs.txt"
    shell:
        """
        python {params.ashic} split2chrs \
            --chr1 1 --allele1 3 \
            --chr2 4 --allele2 6 \
            {input.read_pair} \
            {output} \
            >{log} 2>&1
        """

rule dip3d_ashic_binning:
    input:
        read_pair_dir = rules.dip3d_ashic_split2chrs.output,
        chr_sizes = rules.dip3d_get_chr_sizes.output.txt
    output:
        directory("analysis_other/dip3d/{asm}/5-ashic/{chr}/binned")
    conda:
        "../env/ashic.yml"
    params:
        ashic_data = config['ashic_data'],
        chr = lambda wc: wc.chr,
        res = config.get('ashic_resolution', 10000)
    log:
        "logs/dip3d/ashic_binning/{asm}.{chr}.ashic_binning.log"
    benchmark:
        "runtimes/dip3d/dip3d_ashic_binning/{asm}.{chr}.dip3d_ashic_binning.txt"
    shell:
        """
        python {params.ashic_data} binning \
            --c1=1 --p1=2 --a1=3 \
            --c2=4 --p2=5 --a2=6 \
            --res={params.res} \
            --chrom={params.chr} \
            --genome {input.chr_sizes} \
            {input.read_pair_dir}/ashic_read_pair \
            {output} \
            >{log} 2>&1
        """

rule dip3d_ashic_pack:
    input:
        read_pair_dir = rules.dip3d_ashic_binning.output
    output:
        "analysis_other/dip3d/{asm}/5-ashic/{chr}/packed"
    conda:
        "../env/ashic.yml"
    params:
        ashic_data = config['ashic_data']
    log:
        "logs/dip3d/ashic_pack/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_ashic_pack/{asm}.{chr}.dip3d_ashic_pack.txt"
    shell:
        """
        python {params.ashic_data} pack \
            {input.read_pair_dir} \
            {output} \
            >{log} 2>&1 && \
        """

rule dip3d_run_ashic:
    input:
        rules.dip3d_ashic_pack.output
    output:
        result = "analysis_other/dip3d/{asm}/5-ashic/{chr}/ashic_result.txt"
    conda:
        "../env/ashic.yml"
    params:
        ashic = config['ashic'],
        res = config.get('ashic_resolution', 10000),
        model = config.get('ashic_model')
    log:
        "logs/dip3d/run_ashic/{asm}.{chr}.run_ashic.log"
    benchmark:
        "runtimes/dip3d/dip3d_run_ashic/{asm}.{chr}.dip3d_run_ashic.txt"
    shell:
        """
        python {params.ashic} \
            -i {input}/ashic_read_pair_{wildcards.chr}_{params.res}.pickle \
            -o $(dirname {output.result}) \
            --model {params.model} \
            >{log} 2>&1
        """
