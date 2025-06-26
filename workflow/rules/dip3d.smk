def get_chr_list_for_asm(wc):
    
    file = f"assembly/qc/phased_verkko/{wc.asm}/sample_sex.txt"

    if open(file).read().strip() == "male":
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X", "Y"]] 
    else:
        chroms = ['chr'+str(i) for i in list(range(1, 23)) + ["X"]] 
    
    return(chroms)


rule all_merge_dip3d:
    input:
        "assembly/qc/phased_verkko/{asm}/sample_sex.txt",
        lambda wc: expand("analysis_other/dip3d/{{asm}}/ashic/{chr}/", chr = get_chr_list_for_asm(wc))
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
# Using either Falign (Dip3d paper) or minimap2 (wf-pore-c)

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
        fq = "analysis_other/dip3d/{asm}/{asm}.POREC.fastq.gz",
        fa = config['ref'],
        bed = config['ref'] + ".repeat_regions.bed"
    output: 
        bam = "analysis_other/dip3d/{asm}/1-falign/porec.fragments.bam"
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

rule dip3d_split_bams:
    input:
        bam = rules.dip3d_map_reads.output.bam
    output:
        "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam"
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

# Step 2: Create porec frag pairs:
# -------------------------------

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

rule dip3d_select_chr_snp_bam:
    input:
        bam_dirs = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam"
    output:
        bam = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.hapcut.bam"
    params:
        dip3d = config['dip3d'],
        l = 100,
        q = 5,
        i = 90,
        cov = 30
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

rule dip3d_extract_hairs:
    input:
        ref = config['ref'],
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.vcf",
        bam = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.hapcut.bam"
    output:
        frag = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.frag"
    conda:
        "../env/hapcut2.yml"
    log:
        "logs/dip3d/extract_hairs/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_extract_hairs/{asm}.{chr}.dip3d_extract_hairs.txt"
    shell:
        """
        extractHAIRS \
            --ont 1 \
            --mmq 5 \
            --ref {input.ref} \
            --VCF {input.vcf} \
            --bam {input.bam} \
            --out {output.frag} \
            >{log} 2>&1
        """

rule dip3d_make_frag_pairs:
    input:
        frag = rules.dip3d_extract_hairs.output.frag
    output:
        pairs = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.frag-pair"
    params:
        dip3d = config['dip3d']
    log:
        "logs/dip3d/make_frag_pairs/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_make_frag_pairs/{asm}.{chr}.dip3d_make_frag_pairs.txt"
    shell:
        """
        {params.dip3d} make-pore-c-frag-pair \
            {input.frag} \
            {output.pairs} \
            >{log} 2>&1
        """

rule dip3d_hapcut_hap_asm:
    input:
        pairs = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.frag-pair",
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.vcf",
    output:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.hapcut.phased.VCF"
    conda:
        "../env/hapcut2.yml"
    log:
        "logs/dip3d/hapcut_hap_asm/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/hapcut_hap_asm/{asm}.{chr}.hapcut_hap_asm.txt"
    shell:
        """
        HAPCUT2 \
            --fragments {input.pairs} \
			--VCF {input.vcf} \
			--output analysis_other/dip3d/{wildcards.asm}/3-snp/{wildcards.chr}/{wildcards.chr}.hapcut \
			--hic 1 \
			> {log} 2>&1
        """

# Step 2: Haplotag monomers
# ------------------------

rule dip3d_whatshap_phase_vcf:
    input:
        ref = config['ref'],
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.vcf",
        hapcut_vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.hapcut.phased.VCF",
        hapcut_bam = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.hapcut.bam"
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
            {input.vcf} {input.vcf} {input.hapcut_vcf} {input.hapcut_bam}\
            >{log} 2>&1
        """

rule dip3d_extract_mvp_block:
    input:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.snp.whatshap.vcf"
    output:
        vcf = "analysis_other/dip3d/{asm}/4-bam-haplotag/{chr}/{chr}.mvpblock.vcf"
    params:
        dip3d = config['dip3d']
    log:
        "logs/dip3d/extract_mvp_block/{asm}_{chr}.log"
    benchmark:
        "runtimes/dip3d/dip3d_extract_mvp_block/{asm}.{chr}.dip3d_extract_mvp_block.txt"
    shell:
        """
        {params.dip3d} extract-mvp-het-snp \
            {input.vcf} \
            {output.vcf} \
            >{log} 2>&1
        """

rule dip3d_bam_haplotag:
    input:
        ref = config['ref'],
        vcf = "analysis_other/dip3d/{asm}/4-bam-haplotag/{chr}/{chr}.mvpblock.vcf",
        bam = "analysis_other/dip3d/{asm}/2-chr-bams/{chr}/{chr}.bam"
    output:
        bam = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/{chr}.haplotag"
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

## Step 5: Haplotype imputation (ASHIC)

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
        frag_pair = "analysis_other/dip3d/{asm}/4-haplotag/{chr}/{chr}.haplotag"
    output:
        directory("analysis_other/dip3d/{asm}/ashic-read-pair/{chr}")
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
            {input.frag_pair} \
            {params.chr} \
            {params.min_mapQ} \
            {output} \
            >{log} 2>&1
        """
    
rule dip3d_ashic_split2chrs:
    input:
        read_pair = rules.dip3d_frag_to_ashic_read_pair.output
    output:
        directory("analysis_other/dip3d/{asm}/ashic/{chr}/")
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
        {params.ashic} split2chrs \
            --chr1 1 --allele 3 \
            --chr2 4 --allele2 6 \
            {input.read_pair} \
            >{log} 2>&1
        """

rule dip3d_ashic_binning:
    input:
        read_pair_dir = rules.dip3d_ashic_split2chrs.output,
        chr_sizes = rules.dip3d_get_chr_sizes.output.txt
    output:
        binned = "analysis_other/dip3d/{asm}/ashic/{chr}/ashic-read-pair.pickle"
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
        {params.ashic_data} binning \
            --c1=1 --p1=2 --a1=3 \
            --c2=4 --p2=5 --a2=6 \
            --res={params.res} \
            --chrom={params.chr} \
            --genome {input.chr_sizes} \
            {input.read_pair_dir}/ashic-read-pair \
            {input.read_pair_dir} \
            >{log} 2>&1
        """

rule dip3d_ashic_pack:
    input:
        read_pair_dir = rules.dip3d_ashic_split2chrs.output
    output:
        packed = "analysis_other/dip3d/{asm}/ashic/{chr}/packed.done"
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
        {params.ashic_data} pack \
            {input.read_pair_dir} \
            {input.read_pair_dir} \
            >{log} 2>&1 && \
        touch {output.packed}
        """

rule dip3d_run_ashic:
    input:
        binned = rules.dip3d_ashic_binning.output.binned,
        packed = rules.dip3d_ashic_pack.output.packed
    output:
        result = "analysis_other/dip3d/{asm}/ashic/{chr}/ashic_result.txt"
    conda:
        "../env/ashic.yml"
    params:
        ashic = config['ashic'],
        res = config.get('ashic_resolution', 10000),
        model = config.get('ashic_model', 'default')
    log:
        "logs/dip3d/run_ashic/{asm}.{chr}.run_ashic.log"
    benchmark:
        "runtimes/dip3d/dip3d_run_ashic/{asm}.{chr}.dip3d_run_ashic.txt"
    shell:
        """
        {params.ashic} \
            -i {input.binned} \
            -o $(dirname {output.result}) \
            --model {params.model} \
            >{log} 2>&1
        """
