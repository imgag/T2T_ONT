rule all_gqc_hg002:
    input:
        expand("analysis_other/GQC/{asm}/assemblybench", asm=["T2T_HG002"]),
        expand("analysis_other/GQC/{asm}/readbench",    asm=["T2T_HG002"]),
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/summary.tsv",
        "analysis_other/GQC/T2T_HG002/gqc_error_classification/hp_length_summary.tsv",


rule all_extended_qc:
    input:
        # GQC results ! T2T08 is not working , removed !
        expand("analysis_other/GQC/{asm}/assemblybench", asm=[s for s in finished_samples if s != "T2T08"]),
        expand("analysis_other/GQC/{asm}/readbench", asm=[s for s in finished_samples if s != "T2T08"]),
        # Nucflag results
        expand("analysis_other/nucflag/{asm}/nucflag_status.bed", asm=finished_samples),
        # Longdust results
        expand("analysis_other/longdust/{asm}/longdust.bed", asm=finished_samples),
        # RepeatMasker results
        expand("analysis_other/repeatmasker/{asm}/rm_summary/{asm}_sequence_summary.csv", asm=finished_samples),
        # Flagger results
        expand("analysis_other/flagger/{asm}/prediction_summary_final.tsv", asm=finished_samples),
        # Assembly Issue annotation
        expand("analysis_other/annotations/{asm}/assembly_issue_annotations.lifted.gff3", asm=finished_samples),
        # CenMap
        #expand("analysis_other/cenmap/cenmap.done", asm = finished_samples),
        # CenMap Single
        #expand("analysis_other/cenmap_single/{asm}/cenmap.done", asm = finished_samples),
        # TrioEval
        "analysis_other/trioeval/trioeval_summary.tsv",
        "analysis_other/trioeval/trioeval_phaseblocks.tsv",
        # Normalised variant VCFs
        expand("results/{sample}/variants/small_variants.norm.vcf.gz", sample=finished_samples),
        # Mendelian inheritance violations
        expand("analysis_other/mendelian_qc/{sample}/summary.tsv", sample=["T2T00", "T2T04"]),
        # Gap-end and contig-end distances to nearest multicopy gene locus
        expand("analysis_other/multicopy_genes/{asm}/gap_to_multicopy_gene_distance.tsv", asm=finished_samples),
        expand("analysis_other/multicopy_genes/{asm}/contig_ends.bed", asm=finished_samples),
        expand("analysis_other/multicopy_genes/{asm}/contig_ends_to_multicopy_gene_distance.tsv", asm=finished_samples),
        # Small contig alignment regions
        "doc/tables/for_plots/small_contig_satellite_class.all_samples.tsv"
    output:
        "analysis_other/extended_qc.done"
    shell:
        """
        touch {output}
        """

rule all_cenmap:
    input:
        # CenMap
        "analysis_other/cenmap/cenmap.done",
        # CenMap Single
        expand("analysis_other/cenmap_single/{asm}/cenmap.done", asm = finished_samples)


rule all_mendelian_qc:
    input:
        expand(
            [
                "analysis_other/mendelian_qc/{sample}/summary.tsv",
                "analysis_other/mendelian_qc/{sample}/sites.tsv.gz",
                "analysis_other/mendelian_qc/{sample}/merr.vcf.gz",
                "analysis_other/mendelian_qc/{sample}/merr.bed.gz",
            ],
            sample=["T2T00", "T2T04"],
        ),


rule gqc_assemblybench:
    input:
        assembly = lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"],
        ref=config["ref_hg002_q100"]
    output:
        directory("analysis_other/GQC/{asm}/assemblybench")
    conda:
        "../env/gqc.yml"
    log:
        "logs/gqc/{asm}_gqc_assemblybench.log"
    threads: 10
    params:
        config = config['gqc']['config'],
        resourcedir = config['gqc']['resourcedir'],
        venv = config['gqc']['venv'],
        fastk_binfolder = config['gqc']['fastk_binfolder'],
    resources:
        mem_mb=80000  # 80 GB RAM
    shell:
        """
        source {params.venv}/bin/activate
        export PATH=$PATH:{params.fastk_binfolder}
        assemblybench\
            --reffasta {input.ref} \
            --queryfasta {input.assembly} \
            --config {params.config} \
            --prefix {output} \
            -t {threads} \
            --assembly {wildcards.asm} \
            --benchmark HG002v1.1 \
            > {log} 2>&1
        """

rule gqc_readbench:
    input:
        ref=config["ref_hg002_q100"],
        bam = "data/mapped/{asm}/{asm}.HQ.asm.bam"
    output:
        directory("analysis_other/GQC/{asm}/readbench")
    conda:
        "../env/gqc.yml"
    log:
        "logs/gqc/{asm}_gqc_readbench.log"
    threads: 1
    params:
        config = config['gqc']['config'],
        resourcedir = config['gqc']['resourcedir'],
        venv = config['gqc']['venv'],
        fastk_binfolder = config['gqc']['fastk_binfolder'],
    shell:
        """
        source {params.venv}/bin/activate
        export PATH=$PATH:{params.fastk_binfolder}
        readbench\
            --reffasta {input.ref} \
            --bam {input.bam} \
            --config {params.config} \
            --prefix {output} \
            --readsetname {wildcards.asm} \
            >{log} 2>&1
        """


# This is a faster alternative to repeatmasker. Annotates repeats but no classification
rule longdust:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        out="analysis_other/longdust/{asm}/longdust.bed"
    log:
        "logs/longdust/{asm}_longdust.log"
    threads: 1
    benchmark:
        "runtimes/longdust/{asm}.longdust.bed"
    params:
        longdust=config['longdust']
    shell:
        """
        {params.longdust} \
            {input.fa} \
            > {output.out} 2>{log}
        """

# Repeatmaster Quick
rule repeatmasker_quick:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased": "phased"})["assembly"]
    output:
        out="analysis_other/repeatmasker/{asm}/assembly.fasta.out",
        gff="analysis_other/repeatmasker/{asm}/assembly.fasta.out.gff",
        tbl="analysis_other/repeatmasker/{asm}/assembly.fasta.tbl"
    conda:
        "../env/repeatmasker.yml"
    log:
        "logs/repeatmasker/{asm}_repeatmasker.log"
    threads: 32
    benchmark:
        "runtimes/repeatmasker/{asm}.repeatmasker.txt"
    resources:
        mem_mb=80000  # 80 GB RAM
    params:
        outdir=lambda wc: f"analysis_other/repeatmasker/{wc.asm}",
        species="human"  # Use built-in RepeatMasker library for human
    shell:
        """
        set -euo pipefail

        LOG=$(realpath {log})
        FA=$(realpath {input.fa})
        OUTDIR=$(realpath {params.outdir})

        # Create output directory
        mkdir -p $OUTDIR

        # Change to output directory
        cd $OUTDIR

        # Copy input fasta
        cp $FA assembly.fasta

        # Run RepeatMasker with species library (much faster than HMM)
        # -species: uses curated species-specific repeat library
        # -xsmall: returns repetitive regions in lowercase
        # -gff: creates GFF output
        # -pa: number of parallel threads
        RepeatMasker \
            -species {params.species} \
            -pa {threads} \
            -xsmall \
            -gff \
            assembly.fasta \
            >> $LOG 2>&1

        # Remove temporary files
        rm -f assembly.fasta
        rm -f assembly.fasta.cat
        rm -f *.ori.out

        # Clean up RepeatMasker temporary directories
        rm -rf RM_*/

        echo "RepeatMasker completed successfully" >> $LOG
        """

rule analyze_repeatmasker:
    input:
        rm_out="analysis_other/repeatmasker/{asm}/assembly.fasta.out"
    output:
        summary="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_sequence_summary.csv",
        bed = "analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucflag.bed"
    conda:
        "../env/R.yml"
    log:
        "logs/analyze_repeatmasker/{asm}_analyze_repeatmasker.log"
    shell:
        """
        # Run the R script to analyze the RepeatMasker output
        Rscript workflow/scripts/16_analyze_repeatmasker.R \
            --input {input.rm_out} \
            --output $(dirname {output.summary})/{wildcards.asm}   \
            > {log} 2>&1
        """

rule nucflag:
    input:
        bam="data/mapped/{asm}/{asm}.HQ.asm.bam",
        rm_bed="analysis_other/repeatmasker/{asm}/rm_summary/{asm}_nucflag.bed"
    output:
        status="analysis_other/nucflag/{asm}/nucflag_status.bed",
        misasm="analysis_other/nucflag/{asm}/nucflag_misasm.bed",
        plots=directory("analysis_other/nucflag/{asm}/plots")
    conda:
        "../env/nucflag.yml"
    log:
        "logs/nucflag/{asm}_nucflag.log"
    threads: 12
    shell:
        """
        nucflag \
            -i {input.bam} \
            -d {output.plots} \
            --overlay_regions {input.rm_bed} \
            --threads {threads} \
            --output_status {output.status} \
            --output_misasm {output.misasm} \
            > {log} 2>&1

        # Replace colons with underscores in plot filenames
        for file in {output.plots}/*; do
            if [[ -f "$file" ]]; then
                newname=$(echo "$file" | tr ':' '_')
                if [[ "$file" != "$newname" ]]; then
                    mv "$file" "$newname"
                fi
            fi
        done
        """

rule asm_index:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
    output:
        "analysis_other/flagger/{asm}/assembly.index.bed"
    log:
        "logs/flagger/{asm}/asm_index.log"
    threads:
        1
    shell:
        """
        samtools faidx {input.fa} && \
        awk 'BEGIN{{OFS="\t"}} {{print $1, 0, $2}}' {input.fa}.fai > {output} 2>{log}
        """

rule create_index_json:
    input:
        "analysis_other/flagger/{asm}/assembly.index.bed"
    output:
        "analysis_other/flagger/{asm}/assembly.index.json"
    log:
        "logs/flagger/{asm}/create_json.log"
    params:
        workdir = workflow.basedir
    threads:
        1
    shell:
        """
        echo "{{" > {output}
        echo \\"whole_genome\\" : \\"/mnt/{input}\\" >> {output}
        echo "}}" >> {output}
        """

rule flagger_create_cov:
    input:
        bam = "data/mapped/{asm}/{asm}.HQ.asm.bam",
        idx = "analysis_other/flagger/{asm}/assembly.index.json"
    output:
        cov = "analysis_other/flagger/{asm}/coverage_file.cov.gz",
    log:
        "logs/flagger/{asm}/create_cov.log"
    threads:
        16
    params:
        bam2cov = config['bam2cov'],
        workdir = workflow.basedir
    shell:
        """
        docker run --rm \
            --user $(id -u):$(id -g) \
            -v$(dirname {params.workdir}):/mnt \
            {params.bam2cov} \
                --bam /mnt/{input.bam} \
                --output /mnt/{output.cov} \
                --annotationJson /mnt/{input.idx} \
                --threads {threads} \
                --baselineAnnotation "whole_genome" \
            >{log} 2>&1
        """

rule flagger:
    input:
        cov = "analysis_other/flagger/{asm}/coverage_file.cov.gz"
    output:
        tsv = "analysis_other/flagger/{asm}/prediction_summary_final.tsv",
        bed = "analysis_other/flagger/{asm}/final_flagger_prediction.bed"
    log:
        "logs/flagger/{asm}/hmm_flager.log"
    threads:
        16
    params:
        hmm_flagger = config['hmm_flagger'],
        alpha_tsv = config['hmm_flagger_alpha_tsv'],
        workdir = workflow.basedir
    shell:
        """
        docker run --rm \
            --user $(id -u):$(id -g) \
            -v$(dirname {params.workdir}):/mnt \
            {params.hmm_flagger} \
                --input /mnt/{input.cov} \
                --outputDir /mnt/$(dirname {output.tsv}) \
                --trackName {wildcards.asm}_flagger_prediction \
                --threads {threads} \
        >{log} 2>&1
        """

bed_paths = {
    'flagger': "analysis_other/flagger/{asm}/final_flagger_prediction.bed",
    'nucflag': "analysis_other/nucflag/{asm}/nucflag_misasm.bed",
    'gaps': "assembly/qc/phased_verkko/{asm}/gap_stats.both.n_regions.bed"
}

rule liftover_bed_paf:
    input:
        bed=lambda wc: bed_paths[wc.tool].format(asm=wc.asm),
        paf="assembly/qc/phased_verkko/{asm}/both.mapped_T2T.paf",
    output:
        lifted_bed="analysis_other/annotations/{asm}/{tool}.lifted.bed",
        unmapped="analysis_other/annotations/{asm}/{tool}.unmapped.bed"
    conda:
        "../env/py_report.yml"
    log:
        "logs/annotation/{asm}/liftover_{tool}.log"
    threads: 1
    shell:
        """
        python3 workflow/scripts/40_liftover_bed_paf.py \
            --bed {input.bed} \
            --paf {input.paf} \
            --output {output.lifted_bed} \
            --unmapped {output.unmapped} \
            --min-mapq 5 \
            --min-len 50000 \
            --max-search 10000 \
            --window 10000 \
            --debug \
            > {log} 2>&1

        # Log statistics
        orig_lines=$(grep -v '^#\|^track' {input.bed} | wc -l)
        lifted_lines=$(grep -v '^#\|^track' {output.lifted_bed} | wc -l)
        unmapped_lines=$(wc -l < {output.unmapped})

        echo "Original regions: $orig_lines" >> {log}
        echo "Lifted regions: $lifted_lines" >> {log}
        echo "Unmapped regions: $unmapped_lines" >> {log}
        echo "Success rate: $(echo "scale=2; $lifted_lines/($lifted_lines+$unmapped_lines)*100" | bc)%" >> {log}
        """

rule gff_from_lifted_flagger_nucflag:
    input:
        flagger="analysis_other/annotations/{asm}/flagger.lifted.bed",
        nucflag="analysis_other/annotations/{asm}/nucflag.lifted.bed",
        gaps="analysis_other/annotations/{asm}/gaps.lifted.bed"
    output:
        gff="analysis_other/annotations/{asm}/assembly_issue_annotations.lifted.gff3"
    conda:
        "../env/py_report.yml"
    log:
        "logs/annotation/{asm}/gff_from_lifted_bed.log"
    shell:
        """
        bash workflow/scripts/35_combine_bed_to_gff3.sh \
            {input.flagger} \
            {input.nucflag} \
            {input.gaps} \
            {output.gff} \
            > {log} 2>&1
        """

rule gff_from_flagger_nucflag:
    input:
        flagger="analysis_other/flagger/{asm}/final_flagger_prediction.bed",
        nucflag="analysis_other/nucflag/{asm}/nucflag_status.bed",
        gaps="assembly/qc/phased_verkko/{asm}/gap_stats.both.n_regions.bed"
    output:
        gff="analysis_other/annotations/{asm}/flagger_nucflag_annotations.gff3"
    log:
        "logs/annotation/{asm}/gff_from_flagger_nucflag.log"
    threads: 1
    shell:
        """
        bash workflow/scripts/35_combine_bed_to_gff3.sh \
            {input.flagger} \
            {input.nucflag} \
            {input.gaps} \
            {output.gff} \
            > {log} 2>&1
        """

rule multicopy_gene_loci:
    """Genes overlapping segmental duplications in CHM13 coordinates = candidate multicopy gene loci.

    Reviewer R3 contiguity metric: identifies multicopy gene loci genome-wide
    by intersecting CHM13 gene annotations with the CHM13 SEDEF segmental
    duplication track (Vollger et al. 2022). Same for all samples.
    """
    input:
        genes=config["gene_annotation_chm13"],
        segdup=config["segdup_chm13"],
    output:
        bed="analysis_other/multicopy_genes/chm13_multicopy_genes.bed",
    log:
        "logs/multicopy_genes/intersect.log",
    conda:
        "../env/hapdiff.yml"
    shell:
        """
        zcat {input.genes} \
            | awk -F'\t' 'BEGIN{{OFS="\t"}} $3=="gene" {{
                name="NA";
                n=split($9,a,";");
                for(i=1;i<=n;i++){{ if(a[i] ~ /^gene_name=/){{ split(a[i],b,"="); name=b[2] }} }}
                print $1, $4-1, $5, name
            }}' \
            | sort -k1,1 -k2,2n \
            | bedtools intersect -a - -b {input.segdup} -u \
            > {output.bed} 2>{log}
        """


rule gap_to_multicopy_gene_distance:
    """Distance from each lifted assembly gap end to the nearest multicopy gene locus.

    Reviewer R3 contiguity metric: how close are unresolved contig ends/gaps
    to multicopy genes, in CHM13 coordinates.
    """
    input:
        gaps="analysis_other/annotations/{asm}/gaps.lifted.bed",
        multicopy="analysis_other/multicopy_genes/chm13_multicopy_genes.bed",
    output:
        tsv="analysis_other/multicopy_genes/{asm}/gap_to_multicopy_gene_distance.tsv",
    log:
        "logs/multicopy_genes/{asm}/closest.log",
    conda:
        "../env/hapdiff.yml"
    shell:
        """
        echo -e "chrom\tgap_start\tgap_end\tgap_type\tgap_size\tstrand\tnearest_gene_chrom\tnearest_gene_start\tnearest_gene_end\tnearest_gene_name\tdistance" \
            > {output.tsv}
        sort -k1,1 -k2,2n {input.gaps} \
            | bedtools closest -d -a - -b {input.multicopy} \
            >> {output.tsv} 2>{log}
        """


rule contig_ends_bed:
    """Extract alignment end-positions of incomplete (non-chromosome-sized) contigs in CHM13 coordinates.

    Reviewer R3 contiguity metric: where do smaller, fragmented contigs terminate
    on the CHM13 reference? PAF target coords are already in CHM13 space so no
    liftover step is needed.  Per contig + chromosome, we emit 1-bp BED entries
    at the leftmost and rightmost alignment boundaries.

    Filters applied:
      - mapq >= 5 and alignment block length >= 10 kb
      - contig length < 40 Mb  (excludes chromosome-spanning T2T contigs; the
        smallest human chromosome, chr21, is ~45 Mb)
      - contig length >= 10 kb
    """
    input:
        paf="assembly/qc/phased_verkko/{asm}/both.mapped_T2T.paf",
    output:
        bed="analysis_other/multicopy_genes/{asm}/contig_ends.bed",
    log:
        "logs/multicopy_genes/{asm}/contig_ends.log",
    conda:
        "../env/hapdiff.yml"
    shell:
        """
        awk -F'\t' '
            $12 >= 5 && $11 >= 10000 && $2 >= 10000 && $2 < 40000000 {{
                key = $1 "\t" $6
                if (!(key in min_ts) || $8 < min_ts[key])  min_ts[key] = $8
                if (!(key in max_te) || $9 > max_te[key])  max_te[key] = $9
                clen[key] = $2
                chrom[key] = $6
                hap[key]   = ($1 ~ /haplotype1/) ? "haplotype1" : \
                             ($1 ~ /haplotype2/) ? "haplotype2" : "unassigned"
                cname[key] = $1
            }}
            END {{
                for (k in min_ts) {{
                    s = min_ts[k]
                    e = max_te[k]
                    chr = chrom[k]
                    print chr, s,   s+1, cname[k], clen[k], hap[k], "left_end"
                    print chr, e-1, e,   cname[k], clen[k], hap[k], "right_end"
                }}
            }}
        ' OFS='\t' {input.paf} \
            | sort -k1,1 -k2,2n \
            > {output.bed} 2>{log}
        echo "Contig ends extracted: $(wc -l < {output.bed})" >> {log}
        """


rule contig_ends_to_multicopy_gene_distance:
    """Distance from each incomplete contig end to the nearest multicopy gene locus.

    Reviewer R3 contiguity metric: how close are the ends of smaller, fragmented
    contigs to multicopy gene regions in CHM13 coordinates.
    """
    input:
        ends="analysis_other/multicopy_genes/{asm}/contig_ends.bed",
        multicopy="analysis_other/multicopy_genes/chm13_multicopy_genes.bed",
    output:
        tsv="analysis_other/multicopy_genes/{asm}/contig_ends_to_multicopy_gene_distance.tsv",
    log:
        "logs/multicopy_genes/{asm}/contig_ends_closest.log",
    conda:
        "../env/hapdiff.yml"
    shell:
        """
        echo -e "chrom\tend_start\tend_end\tcontig_name\tcontig_length\thaplotype\tend_side\tnearest_gene_chrom\tnearest_gene_start\tnearest_gene_end\tnearest_gene_name\tdistance" \
            > {output.tsv}
        bedtools closest -d -a {input.ends} -b {input.multicopy} \
            >> {output.tsv} 2>{log}
        """


# rule liftoff_gff:
#     input:
#         gff="analysis_other/annotations/{asm}/flagger_nucflag_annotations.gff3",
#         ref=config["ref"],
#         asm=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"]
#     output:
#         lifted_gff="analysis_other/annotations/{asm}/flagger_nucflag_annotations.lifted.gff3",
#         unmapped = "analysis_other/annotations/{asm}/lifted_unmapped.txt",
#     conda:
#         "../env/liftoff.yml"
#     log:
#         "logs/annotation/{asm}/liftoff_gff.log"
#     threads: 8
#     shell:
#         """
#         feature_file=$(mktemp)
#         cut -f3 {input.gff} | sort | uniq | grep -v '##' > ${{feature_file}}

#         liftoff \
#             -f ${{feature_file}} \
#             -g {input.gff} \
#             -o {output.lifted_gff} \
#             -dir $(dirname {output.lifted_gff}) \
#             -p {threads} \
#             -u {output.unmapped} \
#             {input.ref} {input.asm} \
#             > {log} 2>&1
#         """

rule small_contig_paf_to_bed:
    """Convert PAF to reference-coordinate BED for small contigs.

    Filters:
      - contig length: 10 kb–10 Mb  (10 kb floor; 10 Mb ceiling excludes
        chromosome-spanning T2T contigs)
      - mapq ≥ 5
      - alignment block ≥ 10 kb
    Output columns: chrom, start, end, contig_name, contig_length_bp, haplotype
    """
    input:
        paf="assembly/qc/phased_verkko/{asm}/both.mapped_T2T.paf",
    output:
        bed="analysis_other/multicopy_genes/{asm}/small_contig_alignments.bed",
    log:
        "logs/multicopy_genes/{asm}/small_contig_paf_to_bed.log",
    shell:
        """
        awk -F'\t' '
            $2 >= 10000 && $2 < 10000000 && $12 >= 5 && ($9-$8) >= 10000 {{
                hap = ($1 ~ /haplotype1/) ? "haplotype1" : \
                      ($1 ~ /haplotype2/) ? "haplotype2" : "unassigned"
                print $6, $8, $9, $1, $2, hap
            }}
        ' OFS='\t' {input.paf} | sort -k1,1 -k2,2n > {output.bed} 2>{log}
        echo "Alignment records: $(wc -l < {output.bed})" >> {log}
        """


rule censat_strip_header:
    """Strip track/browser header from the CHM13 CenSat BED (reference, runs once).

    Keeps only chrom/start/end and sorts for bedtools. The raw file has a UCSC
    track header line that bedtools does not accept.
    """
    input:
        bed=config["censat_chm13"],
    output:
        bed="data/ref/chm13v2.0_censat_v2.0.stripped.bed",
    log:
        "logs/multicopy_genes/censat_strip_header.log",
    shell:
        """
        grep -v '^track\|^browser\|^#' {input.bed} \
            | cut -f1-3 | sort -k1,1 -k2,2n \
            > {output.bed} 2>{log}
        echo "CenSat regions: $(wc -l < {output.bed})" >> {log}
        """


rule small_contig_satellite_class:
    """Classify small contigs by satellite fraction of aligned reference span.

    Reviewer R3 contiguity metric: report the number of non-satellite small contigs
    separately from satellite ones, since the latter are expected artefacts of
    assembling repetitive pericentromeric/centromeric sequence.

    Intersects small-contig reference alignments with the CHM13v2.0 CenSat
    annotation. Contigs where ≥50% of aligned bases overlap satellite sequence
    are labelled 'satellite'; the rest 'non_satellite'.

    Output TSV columns:
      contig_name, contig_length_bp, haplotype, primary_chrom,
      aligned_bp, satellite_bp, satellite_fraction, classification
    """
    input:
        aln_bed="analysis_other/multicopy_genes/{asm}/small_contig_alignments.bed",
        sat_bed="data/ref/chm13v2.0_censat_v2.0.stripped.bed",
    output:
        tsv="analysis_other/multicopy_genes/{asm}/small_contig_satellite_class.tsv",
    log:
        "logs/multicopy_genes/{asm}/small_contig_sat_class.log",
    conda:
        "../env/hapdiff.yml"
    shell:
        """
        bedtools intersect \
            -a {input.aln_bed} \
            -b {input.sat_bed} \
            -wao \
            | awk -F'\t' '
                BEGIN {{ OFS="\t" }}
                {{
                    key = $4
                    if (!(key in clen)) {{
                        clen[key] = $5
                        hap[key]  = $6
                    }}
                    aln_len = $3 - $2
                    sat_bp  = $NF + 0
                    chrom_aln[key][$1] += aln_len
                    total_aln[key]     += aln_len
                    total_sat[key]     += sat_bp
                }}
                END {{
                    print "contig_name", "contig_length_bp", "haplotype", "primary_chrom", \
                          "aligned_bp", "satellite_bp", "satellite_fraction", "classification"
                    for (k in total_aln) {{
                        pchrom = ""; max_bp = 0
                        for (c in chrom_aln[k]) {{
                            if (chrom_aln[k][c] > max_bp) {{ max_bp = chrom_aln[k][c]; pchrom = c }}
                        }}
                        frac  = (total_aln[k] > 0) ? total_sat[k] / total_aln[k] : 0
                        class = (frac >= 0.5) ? "satellite" : "non_satellite"
                        print k, clen[k], hap[k], pchrom, \
                              total_aln[k], total_sat[k], frac, class
                    }}
                }}
            ' OFS='\t' > {output.tsv} 2>{log}
        echo "Small contigs classified: $(tail -n+2 {output.tsv} | wc -l)" >> {log}
        """


rule collect_small_contig_satellite_class:
    """Combine per-sample small-contig satellite classifications into one TSV."""
    input:
        expand(
            "analysis_other/multicopy_genes/{asm}/small_contig_satellite_class.tsv",
            asm=finished_samples,
        ),
    output:
        tsv="doc/tables/for_plots/small_contig_satellite_class.all_samples.tsv",
    log:
        "logs/multicopy_genes/collect_small_contig_sat_class.log",
    shell:
        """
        (
            head -1 {input[0]} | awk 'BEGIN{{OFS="\\t"}} {{print "sample", $0}}'
            for f in {input}; do
                sample=$(basename $(dirname "$f"))
                tail -n +2 "$f" | awk -v s="$sample" 'BEGIN{{OFS="\\t"}} {{print s, $0}}'
            done
        ) > {output.tsv} 2>{log}
        """


rule link_cenmap_input:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"],
        hq="assembly/input/{asm}/{asm}.HQ_herro.50x.fastq.gz",
        mod=lambda wc: find_input_datasets(SimpleNamespace(dataset=wc.asm, type="UL"))["files"][0]
    output:
        fa_link="analysis_other/cenmap/data/asm/{asm}/assembly.fasta",
        hq_link="analysis_other/cenmap/data/hifi/{asm}/{asm}.HQ_herro.50x.fastq.gz",
        hq_fofn="analysis_other/cenmap/data/hifi/{asm}.fofn",
        mod_link="analysis_other/cenmap/data/ont/{asm}/{asm}.unmapped.mod.bam"
    shell:
        """
        mkdir -p $(dirname {output.fa_link})
        ln -sf $(realpath {input.fa}) {output.fa_link}

        mkdir -p $(dirname {output.hq_link})
        ln -sf $(realpath {input.hq}) {output.hq_link}

        echo $(realpath {input.hq}) > {output.hq_fofn}

        mkdir -p $(dirname {output.mod_link})
        ln -sf $(realpath {input.mod}) {output.mod_link}
        """

rule link_cenmap_single_input:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased" : "phased"})["assembly"],
        hq="assembly/input/{asm}/{asm}.HQ_herro.50x.fastq.gz",
        mod=lambda wc: find_input_datasets(SimpleNamespace(dataset=wc.asm, type="UL"))["files"][0]
    output:
        fa_link="analysis_other/cenmap_single/{asm}/data/asm/{asm}/assembly.fasta",
        hq_link="analysis_other/cenmap_single/{asm}/data/hifi/{asm}/{asm}.HQ_herro.50x.fastq.gz",
        hq_fofn="analysis_other/cenmap_single/{asm}/data/hifi/{asm}.fofn",
        mod_link="analysis_other/cenmap_single/{asm}/data/ont/{asm}/{asm}.unmapped.mod.bam"
    shell:
        """
        mkdir -p $(dirname {output.fa_link})
        ln -sf $(realpath {input.fa}) {output.fa_link}

        mkdir -p $(dirname {output.hq_link})
        ln -sf $(realpath {input.hq}) {output.hq_link}

        echo $(realpath {input.hq}) > {output.hq_fofn}

        mkdir -p $(dirname {output.mod_link})
        ln -sf $(realpath {input.mod}) {output.mod_link}
        """

# Cenmap
rule cenmap:
    input:
        fa_link=expand("analysis_other/cenmap/data/asm/{asm}/assembly.fasta", asm = finished_samples),
        hq_link=expand("analysis_other/cenmap/data/hifi/{asm}/{asm}.HQ_herro.50x.fastq.gz", asm = finished_samples),
        hq_fofn=expand("analysis_other/cenmap/data/hifi/{asm}.fofn", asm = finished_samples),
        mod_link=expand("analysis_other/cenmap/data/ont/{asm}/{asm}.unmapped.mod.bam", asm = finished_samples)
    output:
        done = "analysis_other/cenmap/cenmap.done"
    conda:
        "../env/cenmap.yml"
    log:
        "logs/cenmap.log"
    benchmark:
        "runtimes/cenmap/cenmap.txt"
    threads:
        96
    shell:
        """
        LOG=$(realpath {log})
        WD=$(dirname {output.done})
        pushd $WD >{log}

        snakemake \
            --configfile config.yaml \
            --cores {threads} \
            >$LOG 2>&1
        touch cenmap.done
        """

rule cenmap_single_config:
    output:
        cfg          = "analysis_other/cenmap_single/{asm}/config.yaml",
        nucflag      = "analysis_other/cenmap_single/{asm}/config/nucflag.toml",
        rm_colors    = "analysis_other/cenmap_single/{asm}/config/repeatmasker_colors.json",
        rm_sat       = "analysis_other/cenmap_single/{asm}/config/repeatmasker_sat_annot_colors.json",
        stv_colors   = "analysis_other/cenmap_single/{asm}/config/stv_annot_colors.tsv",
    params:
        config   = config["cenmap_single_configfile"],
        cfg_src  = "bin/cenmap/config",
    shell:
        """
        echo 'samples:' > {output.cfg}
        echo '  - {wildcards.asm}' >> {output.cfg}
        echo ' ' >> {output.cfg}
        cat {params.config} >> {output.cfg}

        mkdir -p $(dirname {output.nucflag})
        ln -sf $(realpath {params.cfg_src}/nucflag.toml)                       {output.nucflag}
        ln -sf $(realpath {params.cfg_src}/repeatmasker_colors.json)            {output.rm_colors}
        ln -sf $(realpath {params.cfg_src}/repeatmasker_sat_annot_colors.json)  {output.rm_sat}
        ln -sf $(realpath {params.cfg_src}/stv_annot_colors.tsv)                {output.stv_colors}
        """

# Run cenmap isolated separate on every sample
rule cenmap_single:
    input:
        fa_link="analysis_other/cenmap_single/{asm}/data/asm/{asm}/assembly.fasta",
        hq_link="analysis_other/cenmap_single/{asm}/data/hifi/{asm}/{asm}.HQ_herro.50x.fastq.gz",
        hq_fofn="analysis_other/cenmap_single/{asm}/data/hifi/{asm}.fofn",
        mod_link="analysis_other/cenmap_single/{asm}/data/ont/{asm}/{asm}.unmapped.mod.bam",
        config=ancient("analysis_other/cenmap_single/{asm}/config.yaml")
    output:
        done = "analysis_other/cenmap_single/{asm}/cenmap.done",
        bed_c = "analysis_other/cenmap_single/{asm}/results/final/bed/{asm}_complete_cens.bed",
        bed_cc = "analysis_other/cenmap_single/{asm}/results/final/bed/{asm}_complete_correct_cens.bed",


    conda:
        "../env/cenmap.yml"
    log:
        "logs/cenmap_single/{asm}.log"
    benchmark:
        "runtimes/cenmap_single/{asm}/cenmap.txt"
    threads:
        24
    params:
        config = config["cenmap_single_configfile"]
    shell:
        """
        cenmap \
            --config {input.config} \
            --processes {threads} \
            --output-dir "analysis_other/cenmap_single/{wildcards.asm}" \
            >{log} 2>&1

        touch {output}   
        """


## Trio comparison - evaluate each haplotype separately
rule all_trioeval:
    input:
        "analysis_other/trioeval/trioeval_summary.tsv",
        "analysis_other/trioeval/trioeval_phaseblocks.tsv",


## Trio comparison - evaluate each haplotype separately
rule yak_trioeval_per_haplotype:
    input:
        asm = lambda wc: (
            f"assembly/output/final/{wc.asm.replace('_pofo', '')}/assembly.hap{wc.hp}.pansn.fasta"
            if wc.asm.endswith("_pofo")
            else (f"assembly/output/hifiasm/{wc.asm}/assembly.haplotype{wc.hp}.fasta"
                  if "hifiasm" in wc.asm
                  else f"assembly/output/verkko/{wc.asm}/assembly.haplotype{wc.hp}.fasta")
        ),
        yak_pat = lambda wc: "assembly/input/{}/paternal.yak".format(
            asm[wc.asm.replace("_pofo", "")]["dataset"] if wc.asm.endswith("_pofo") else asm[wc.asm]["dataset"]
        ),
        yak_mat = lambda wc: "assembly/input/{}/maternal.yak".format(
            asm[wc.asm.replace("_pofo", "")]["dataset"] if wc.asm.endswith("_pofo") else asm[wc.asm]["dataset"]
        ),
    output:
        "analysis_other/trioeval/{asm}.haplotype{hp}.trio_eval.txt",
    params:
        yak = config['yak'],
    log:
        "logs/trioeval/{asm}/haplotype{hp}.log",
    shell:
        """
        {params.yak} trioeval \
            {input.yak_pat} {input.yak_mat} {input.asm} \
            > {output} 2>{log}
        """


rule parse_trioeval:
    input:
        expand(
            "analysis_other/trioeval/{asm}.haplotype{hp}.trio_eval.txt",
            asm=asm_method_comparison,
            hp=["1", "2"],
        ),
    output:
        summary   = "analysis_other/trioeval/trioeval_summary.tsv",
        per_contig = "analysis_other/trioeval/trioeval_per_contig.tsv",
    params:
        min_contig_length = 500000,
    conda:
        "../env/py_report.yml"
    log:
        "logs/trioeval/parse_trioeval.log",
    shell:
        """
        python3 workflow/scripts/45_parse_trioeval.py \
            --input {input} \
            --summary {output.summary} \
            --per-contig {output.per_contig} \
            --min-length {params.min_contig_length} \
            > {log} 2>&1
        """


rule parse_trioeval_phaseblocks:
    input:
        trioeval = expand("analysis_other/trioeval/{{asm}}.haplotype{hp}.trio_eval.txt", hp=["1", "2"]),
        paf = lambda wc: (
            expand("assembly/qc/phased_verkko/{sample}/haplotype{hp}.mapped_T2T.paf",
                   sample=wc.asm.replace("_pofo", ""), hp=["1", "2"])
            if wc.asm.endswith("_pofo")
            else (expand("assembly/qc/phased_hifiasm/{asm}/haplotype{hp}.mapped_T2T.paf", asm=wc.asm, hp=["1", "2"])
                  if "hifiasm" in wc.asm
                  else expand("assembly/qc/phased_verkko/{asm}/haplotype{hp}.mapped_T2T.paf", asm=wc.asm, hp=["1", "2"]))
        ),
    output:
        phaseblocks = "analysis_other/trioeval/{asm}.phaseblocks.tsv",
        switches    = "analysis_other/trioeval/{asm}.switches.tsv",
    params:
        min_kmer_support = 5,
    conda:
        "../env/py_report.yml"
    log:
        "logs/trioeval/{asm}/parse_phaseblocks.log",
    shell:
        """
        python3 workflow/scripts/46_parse_trioeval_phaseblocks.py \
            --input {input.trioeval} \
            --paf {input.paf} \
            --phaseblocks {output.phaseblocks} \
            --switches {output.switches} \
            --min-kmer-support {params.min_kmer_support} \
            > {log} 2>&1
        """


rule collect_trioeval_phaseblocks:
    input:
        phaseblocks = expand("analysis_other/trioeval/{asm}.phaseblocks.tsv", asm=asm_method_comparison),
        switches    = expand("analysis_other/trioeval/{asm}.switches.tsv",    asm=asm_method_comparison),
    output:
        phaseblocks = "analysis_other/trioeval/trioeval_phaseblocks.tsv",
        switches    = "analysis_other/trioeval/trioeval_switches.tsv",
    shell:
        """
        head -1 {input.phaseblocks[0]} > {output.phaseblocks}
        tail -n +2 -q {input.phaseblocks} >> {output.phaseblocks}

        head -1 {input.switches[0]} > {output.switches}
        tail -n +2 -q {input.switches} >> {output.switches}
        """


## ============================================================
## Mendelian inheritance violation (MIV) analysis
##
## For each trio (T2T00, T2T04) compare the child's phased
## assembly VCF against both parental phased VCFs.
##
## Site selection: biallelic het (GT = 0|1 or 1|0) with AD = 1,1.
## The AD = 1,1 filter is the formal implementation of the reviewer's
## "exactly two contigs" requirement: in dipcall output AD encodes
## the number of assembly contigs supporting each allele; any site
## with AD > 1,1 has ≥3 contigs and is excluded as a potential
## collapsed-repeat / misalignment artefact.
## ============================================================

TRIO_FAMILIES = {
    "T2T00": ("T2T00_1", "T2T00_2"),
    "T2T04": ("T2T04_1", "T2T04_2"),
}

rule mendelian_merge:
    """Merge trio VCFs and apply the two-contig filter.

    Child sites are restricted to biallelic hets (GT = 0|1 or 1|0) with
    AD = 1,1 — exactly one assembly contig per allele.
    """
    input:
        child   = "results/{sample}/variants/small_variants.norm.vcf.gz",
        parent1 = lambda wc: f"results/{TRIO_FAMILIES[wc.sample][0]}/variants/small_variants.norm.vcf.gz",
        parent2 = lambda wc: f"results/{TRIO_FAMILIES[wc.sample][1]}/variants/small_variants.norm.vcf.gz",
    output:
        vcf = temp("analysis_other/mendelian_qc/{sample}/trio_filtered.vcf.gz"),
    log:
        "logs/mendelian_qc/{sample}.merge.log",
    wildcard_constraints:
        sample = "|".join(TRIO_FAMILIES.keys()),
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools merge --missing-to-ref \
                {input.child} {input.parent1} {input.parent2} 2>{log} |
            bcftools view \
                -i '(GT[0]="0|1" || GT[0]="1|0") && AD[0:0]=1 && AD[0:1]=1' \
                -Oz -o {output.vcf} 2>>{log}
        tabix {output.vcf}
        """


rule mendelian_count:
    """Run bcftools +mendelian2 in count mode (-m c) to produce the text stats summary."""
    input:
        vcf = "analysis_other/mendelian_qc/{sample}/trio_filtered.vcf.gz",
    output:
        stats   = "analysis_other/mendelian_qc/{sample}/stats.txt",
        summary = "analysis_other/mendelian_qc/{sample}/summary.tsv",
    params:
        child  = lambda wc: wc.sample,
        father = lambda wc: TRIO_FAMILIES[wc.sample][0],
        mother = lambda wc: TRIO_FAMILIES[wc.sample][1],
    log:
        "logs/mendelian_qc/{sample}.count.log",
    wildcard_constraints:
        sample = "|".join(TRIO_FAMILIES.keys()),
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools +mendelian2 {input.vcf} \
            -p {params.child},{params.father},{params.mother} \
            -m c > {output.stats} 2>{log}

        awk -v s={wildcards.sample} '
            /^nmerr/    {{ nmerr=$2 }}
            /^ngood/    {{ ngood=$2 }}
            /^nmissing/ {{ nmissing=$2 }}
            END {{
                total = nmerr + ngood
                print "sample\\tnmerr\\tngood\\tnmissing\\tmerr_rate"
                printf "%s\\t%d\\t%d\\t%d\\t%.6f\\n", s, nmerr, ngood, nmissing, \
                    (total > 0 ? nmerr/total : 0)
            }}
        ' {output.stats} > {output.summary}
        """


rule mendelian_annotate:
    """Annotate all evaluated sites with MERR count and extract to TSV for R plotting."""
    input:
        vcf = "analysis_other/mendelian_qc/{sample}/trio_filtered.vcf.gz",
    output:
        sites = "analysis_other/mendelian_qc/{sample}/sites.tsv.gz",
    params:
        child  = lambda wc: wc.sample,
        father = lambda wc: TRIO_FAMILIES[wc.sample][0],
        mother = lambda wc: TRIO_FAMILIES[wc.sample][1],
    log:
        "logs/mendelian_qc/{sample}.annotate.log",
    wildcard_constraints:
        sample = "|".join(TRIO_FAMILIES.keys()),
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools +mendelian2 {input.vcf} \
            -p {params.child},{params.father},{params.mother} \
            -m a 2>{log} \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MERR\n' \
        | bgzip > {output.sites}
        """


rule mendelian_error_vcf:
    """Extract error sites (MERR > 0) as an indexed VCF for genome-browser tracks."""
    input:
        vcf = "analysis_other/mendelian_qc/{sample}/trio_filtered.vcf.gz",
    output:
        merr_vcf = "analysis_other/mendelian_qc/{sample}/merr.vcf.gz",
    params:
        child  = lambda wc: wc.sample,
        father = lambda wc: TRIO_FAMILIES[wc.sample][0],
        mother = lambda wc: TRIO_FAMILIES[wc.sample][1],
    log:
        "logs/mendelian_qc/{sample}.error_vcf.log",
    wildcard_constraints:
        sample = "|".join(TRIO_FAMILIES.keys()),
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools +mendelian2 {input.vcf} \
            -p {params.child},{params.father},{params.mother} \
            -m ae -Oz -o {output.merr_vcf} 2>{log}
        tabix {output.merr_vcf}
        """


rule mendelian_error_bed:
    """Convert MIV sites to a tabix-indexed BED for IGV / genome-browser tracks.

    Columns: chrom, start (0-based), end, name (REF>ALT), score (MERR count).
    Load merr.bed.gz + merr.bed.gz.tbi directly in IGV alongside other error tracks.
    """
    input:
        vcf = "analysis_other/mendelian_qc/{sample}/merr.vcf.gz",
    output:
        bed = "analysis_other/mendelian_qc/{sample}/merr.bed.gz",
    log:
        "logs/mendelian_qc/{sample}.error_bed.log",
    wildcard_constraints:
        sample = "|".join(TRIO_FAMILIES.keys()),
    conda:
        "../env/bcftools.yml"
    shell:
        """
        bcftools query -f '%CHROM\t%POS0\t%END\t%REF>%ALT\t%MERR\n' {input.vcf} \
            | bgzip > {output.bed} 2>{log}
        tabix -p bed {output.bed}
        """


# ── Reviewer QC tables ────────────────────────────────────────────────────────

REVIEWER_SAMPLES = [
    "T2T00", "T2T00_1", "T2T00_2",
    "T2T01", "T2T02", "T2T03",
    "T2T04", "T2T04_1", "T2T04_2",
    "T2T05", "T2T06", "T2T07", "T2T08", "T2T09", "T2T10",
    "T2T11", "T2T12", "T2T13",
    "T2T15", "T2T16", "T2T17", "T2T18", "T2T19",
]


def _ul_processed_summary_files(wildcards):
    d = datasets.get(wildcards.sample, {})
    ul = d.get("UL", [])
    if not isinstance(ul, list):
        ul = [ul]
    return [
        f"data/bamstats/basecalled/sup/{Path(p).name}/{Path(p).name}.sup.sequencing_summary.processed.txt"
        for p in ul
    ]


rule raw_ul_stats:
    """Aggregate raw UL read stats (yield, N50, frac >100 kb) for one sample."""
    input:
        _ul_processed_summary_files,
    output:
        "data/qc/raw_ul_stats/{sample}.raw_ul_stats.tsv",
    log:
        "logs/raw_ul_stats/{sample}.log",
    wildcard_constraints:
        sample="|".join(REVIEWER_SAMPLES),
    conda:
        "../env/python.yml"
    script:
        "../scripts/65a_aggregate_raw_ul_stats.py"


rule reviewer_qc_tables:
    """Generate doc/tables/reviewer_sample_qc.{md,docx} and reviewer_sample_qc_wide.{md,docx}."""
    input:
        raw_stats=expand("data/qc/raw_ul_stats/{sample}.raw_ul_stats.tsv", sample=REVIEWER_SAMPLES),
        datasets="data/datasets.yml",
        asm_qc="doc/tables/qc_samples_full_wide.tsv",
    output:
        t1_md="doc/tables/reviewer_sample_qc.md",
        t1_docx="doc/tables/reviewer_sample_qc.docx",
        t2_md="doc/tables/reviewer_sample_qc_wide.md",
        t2_docx="doc/tables/reviewer_sample_qc_wide.docx",
    log:
        "logs/reviewer_qc_tables.log",
    conda:
        "../env/python.yml"
    script:
        "../scripts/65_reviewer_sample_qc_table.py"
