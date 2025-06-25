chrom_list = ['chr'+str(i) for i in range(1, 23)] + ["X", "Y"]

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
    shell:
        """
        {params.falign} \
            -repeat_bed {input.bed} \
            -num_threads {threads} \
            -outfmt frag-bam \
            -out {output.bam} \
            {input.fq} \
            {input.fa} \
            >{log} 2>&1
        """

rule dip3d_split_bams:
    input:
        bam = rules.dip3d_map_reads.output.bam
    output:
        directory("analysis_other/dip3d/{asm}/2-chr-bams/{chr}")
    params:
        dip3d = config['dip3d'],
        chrom_list = chrom_list
    threads: 16
    log:
        "logs/dip3d/split_bams/{asm}_{chr}.log"
    shell:
        """
        {params.dip3d} split-bam \
            {input.bam} \
            $(dirname {output}) \
            {params.chrom_list} \
            >{log} 2>&1
        """

# Step 2: Create porec frag pairs:
# -------------------------------

rule dip3d_split_vcf:
    input:
        vcf = "assembly/variants/{asm}/phased_verkko/small_variants.dip.vcf.gz"
    output:
        vcf = "analysis_other/dip3d/{asm}/3-snp/{chr}/snp.vcf"
    params:
        dip3d = config['dip3d']
    log:
        "logs/dip3d/split_vcf/{asm}_{chr}.log"
    shell:
        """
        {params.dip3d} split-vcf \
            -v -p \
            {input.vcf} \
            $(dirname {output.vcf}) \
            >{log} 2>&1
        """

rule dip3d_select_chr_snp_bam:
    input:
        bam_dir = rules.dip3d_split_bams.output
    output:
        bam = "analysis_other/dip3d/{asm}/3-snp/{chr}/sorted.snp.bam"
    params:
        dip3d = config['dip3d'],
        l = 100,
        q = 5,
        i = 90,
        cov = 30
    threads: 16
    log:
        "logs/dip3d/select_chr_snp_bam/{asm}_{chr}.log"
    shell:
        """
        {params.dip3d} select-chr-snp-bam \
            -l {params.l} \
            -q {params.q} \
            -i {params.i} \
            {input.bam_dir} \
            {params.cov} \
            {output.bam} \
            {wildcards.chr} \
            >{log} 2>&1
        """

rule dip3d_extract_hairs:
    input:
        ref = config['ref'],
        vcf = "assembly/variants/{asm}/phased_verkko/small_variants.dip.vcf.gz",
        bam = rules.dip3d_select_chr_snp_bam.output.bam
    output:
        frag = "analysis_other/dip3d/{asm}/3-snp/{chr}/{chr}.frag"
    conda:
        "../env/hapcut2.yml"
    log:
        "logs/dip3d/extract_hairs/{asm}_{chr}.log"
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
    shell:
        """
        {params.dip3d} make-pore-c-frag-pair \
            {input.frag} \
            {output.pairs} \
            >{log} 2>&1
        """

# Step 2: Haplotag monomers
# ------------------------

rule dip3d_extract_mvp_block:
    input:
        vcf = "assembly/variants/{asm}/phased_verkko/small_variants.dip.vcf.gz"
    output:
        vcf = "analysis_other/dip3d/{asm}/{chr}/{chr}.mvpblock.vcf"
    params:
        dip3d = config['dip3d']
    log:
        "logs/dip3d/extract_mvp_block/{asm}_{chr}.log"
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
        vcf = "analysis_other/dip3d/{asm}/mvpblock.vcf",
        bam = "analysis_other/dip3d/{asm}/porec.fragments.bam"
    output:
        bam = "analysis_other/dip3d/{asm}/porec.fragments.haplotagged.bam"
    params:
        dip3d = config['dip3d']
    threads: 8
    log:
        "logs/dip3d/bam_haplotag/{asm}.log"
    shell:
        """
        {params.dip3d} tag-bam \
            $(dirname {output.bam}) \
            {wildcards.asm} \
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
    shell:
        """
        {params.dip3d} extract-ashic-chr-size \
            {input.ref} \
            {output.txt} \
            >{log} 2>&1
        """

rule dip3d_frag_to_ashic_read_pair:
    input:
        frag_pair = rules.dip3d_make_frag_pairs.output.pairs
    output:
        directory("analysis_other/dip3d/{asm}/ashic-read-pair/{chr}")
    params:
        dip3d = config['dip3d'],
        min_mapQ = 5,
        chr = lambda wc: wc.chr
    log:
        "logs/dip3d/frag_to_ashic_read_pair/{asm}_{chr}.log"
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
        "logs/dip3d/ashic_split2chrs/{asm}_{chr}.log"
    shell:
        """
        {params.ashic} split2chrs \
            --chr1 1 --allele 3 \
            --chr2 4 --allele2 6 \
            {input.read_pair} \
            >{log} 2>&1
        """