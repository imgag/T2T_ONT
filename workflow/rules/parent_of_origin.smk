# Parent-of-origin assignment using PatMat
#
# Strategy: our Verkko + Pore-C assembly is already chromosome-scale phased,
# so the assembly-derived phased VCF (dipcall) serves as both the per-variant
# phased input (-v) AND the chromosome-scale phase reference (-stv). PatMat
# then uses methylation at imprinted DMRs in the haplotagged BAM to determine
# which assembly haplotype is maternal and which is paternal.
#
# Targets
#   analysis_other/parent_of_origin/{sample}/{sample}_PofO_Scores.tsv
#   assembly/output/final/{sample}/assembly.maternal.fasta
#   assembly/output/final/{sample}/assembly.paternal.fasta
#   assembly/output/final/{sample}/assembly.unassigned.fasta
#   assembly/output/final/{sample}/assembly.maternal.pansn.fasta    ← PanSN-spec headers
#   assembly/output/final/{sample}/assembly.paternal.pansn.fasta
#   assembly/output/final/{sample}/assembly.unassigned.pansn.fasta
#   assembly/output/final/{sample}/assembly.combined.pansn.fasta    ← all contigs combined
#   assembly/output/final/{sample}/contig_pofo_assignment.tsv       ← per-contig origin table

rule all_parent_of_origin:
    input:
        expand(
            "assembly/output/final/{sample}/assembly.{origin}.pansn.fasta",
            sample = finished_samples,
            origin = ["hap1", "hap2", "combined"],
        )


rule patmat_filter_vcf:
    """Filter VCF to diploid chromosomes only before patmat.

    Males: autosomes only (chrX/Y are hemizygous — no phased hets, patmat cannot use them).
    Females: autosomes + chrX (diploid); chrY is absent from the assembly.
    """
    input:
        vcf = "assembly/variants/{sample}/phased_verkko/small_variants.vcf.gz",
        sex = "assembly/qc/phased_verkko/{sample}/sample_sex.txt",
    output:
        vcf = temp("analysis_other/parent_of_origin/{sample}/small_variants.diploid.vcf.gz"),
    log:
        "logs/parent_of_origin/{sample}.patmat_filter_vcf.log",
    conda:
        "../env/bcftools.yml"
    shell:
        """
        AUTOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,\
chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
        SEX=$(cat {input.sex})
        if [ "$SEX" = "female" ]; then
            CHROMS="${{AUTOSOMES}},chrX"
        else
            CHROMS="${{AUTOSOMES}}"
        fi
        bcftools view -r $CHROMS {input.vcf} -Oz -o {output.vcf} 2>{log}
        tabix {output.vcf}
        """


rule patmat:
    input:
        vcf = "analysis_other/parent_of_origin/{sample}/small_variants.diploid.vcf.gz",
        bam = "data/mapped/{sample}/{sample}.UL.100G.haplotagged.ref.bam",
        ref = config["ref"],
        dmr = config["patmat_dmr"],
    output:
        scores  = "analysis_other/parent_of_origin/{sample}/{sample}_PofO_Scores.tsv",
        vcf_out = "analysis_other/parent_of_origin/{sample}/{sample}_PofO_Assignment.vcf.gz",
        cram    = "analysis_other/parent_of_origin/{sample}/{sample}_PofO_Tagged.cram",
    params:
        outprefix = "analysis_other/parent_of_origin/{sample}/{sample}",
        vcf_raw   = "analysis_other/parent_of_origin/{sample}/{sample}_PofO_Assigned.vcf",
    log:
        "logs/parent_of_origin/{sample}.patmat.log",
    benchmark:
        "runtimes/patmat/{sample}.txt",
    threads: 8
    conda:
        "/mnt/storage5/users/ahgrosc1/environments/envs/patmap"
    shell:
        # The assembly-phased VCF is chromosome-scale by construction (Verkko + PoreC),
        # so it is passed as both -v (per-variant phasing) and -stv (chromosome-scale
        # reference) — no strand-seq data required.
        # PatMat writes _PofO_Assigned.vcf (uncompressed); bgzip+tabix it to produce
        # the declared _PofO_Assignment.vcf.gz output.
        """
        patmat \
            -v   {input.vcf} \
            -stv {input.vcf} \
            -b   {input.bam} \
            -ref {input.ref} \
            -kd  {input.dmr} \
            -o   {params.outprefix} \
            -p   {threads} \
            >>{log} 2>&1
        bgzip -@ {threads} {params.vcf_raw} >>{log} 2>&1
        mv {params.vcf_raw}.gz {output.vcf_out}
        tabix {output.vcf_out} >>{log} 2>&1
        """


rule assign_haplotype_pofo:
    input:
        scores = "analysis_other/parent_of_origin/{sample}/{sample}_PofO_Scores.tsv",
        paf    = "assembly/qc/phased_verkko/{sample}/both.mapped_T2T.paf",
        hap1   = "assembly/output/verkko/{sample}/assembly.haplotype1.fasta",
        hap2   = "assembly/output/verkko/{sample}/assembly.haplotype2.fasta",
        both   = "assembly/output/verkko/{sample}/assembly.fasta",
    output:
        hap1     = "assembly/output/final/{sample}/assembly.hap1.pansn.fasta",
        hap2     = "assembly/output/final/{sample}/assembly.hap2.pansn.fasta",
        combined = "assembly/output/final/{sample}/assembly.combined.pansn.fasta",
        table    = "assembly/output/final/{sample}/contig_pofo_assignment.tsv",
    log:
        "logs/parent_of_origin/{sample}.assign_haplotype_pofo.log",
    shell:
        """
        python workflow/scripts/62_assign_haplotype_pofo.py \
            --scores    {input.scores} \
            --paf       {input.paf} \
            --hap1      {input.hap1} \
            --hap2      {input.hap2} \
            --both      {input.both} \
            --out-hap1  {output.hap1} \
            --out-hap2  {output.hap2} \
            --combined  {output.combined} \
            --table     {output.table} \
            --sample    {wildcards.sample} \
            >{log} 2>&1
        """


# PanSN-spec header format (applied directly by 62_assign_haplotype_pofo.py):
#   {sample}#1#{contig} — maternal (PofO assigned) or haplotype1 (unassigned)
#   {sample}#2#{contig} — paternal (PofO assigned) or haplotype2 (unassigned)
# Example: T2T00#1#haplotype1-0000001

