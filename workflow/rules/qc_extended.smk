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
        # CenMap
        expand("analysis_other/cenmap/{asm}/cenmap.done", asm = finished_samples)

rule all_annotation:
    input:
        #expand("analysis_other/annotations/{asm}/flagger_nucflag_annotations.lifted.gff3", asm=finished_samples),
        expand("analysis_other/annotations/{asm}/assembly_issue_annotations.lifted.gff3", asm=finished_samples)

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

# Alternative: RepeatMasker with custom Dfam HMM library (slower but more comprehensive)
rule repeatmasker_hmm:
    input:
        fa=lambda wc: get_assembly_output({**wc, "tool": "verkko", "hp": "both", "isphased": "phased"})["assembly"]
    output:
        out="analysis_other/repeatmasker_hmm/{asm}/assembly.fasta.out",
        gff="analysis_other/repeatmasker_hmm/{asm}/assembly.fasta.out.gff",
        tbl="analysis_other/repeatmasker_hmm/{asm}/assembly.fasta.tbl"
    conda:
        "../env/repeatmasker.yml"
    log:
        "logs/repeatmasker/{asm}_repeatmasker_hmm.log"
    threads: 32
    benchmark:
        "runtimes/repeatmasker/{asm}.repeatmasker_hmm.txt"
    resources:
        mem_mb=120000  # 120 GB RAM
    params:
        dfam_lib=config.get('dfam_hmm'),
        outdir=lambda wc: f"analysis_other/repeatmasker_hmm/{wc.asm}",
    shell:
        """
        set -euo pipefail

        LOG=$(realpath {log})
        FA=$(realpath {input.fa})
        LIB=$(realpath {params.dfam_lib})
        OUTDIR=$(realpath {params.outdir})

        # Create output directory
        mkdir -p $OUTDIR

        # Change to output directory - RepeatMasker works best in the directory with the input
        cd $OUTDIR

        # Copy input fasta
        cp $FA assembly.fasta

        # Copy HMM library (not symlink - RepeatMasker may need to modify/index it)
        cp $LIB dfam.hmm

        # Run RepeatMasker with HMM library
        # Using hmmer engine (auto-selects nhmmer for DNA)
        # -q: quick search (5-10% less sensitive, but 3-4x faster)
        # -no_is: skips bacterial insertion element check
        # -xsmall: returns repetitive regions in lowercase
        # -gff: creates GFF output
        RepeatMasker \
            -lib dfam.hmm \
            -engine hmmer \
            -pa {threads} \
            -gff \
            -q \
            -no_is \
            -xsmall \
            assembly.fasta \
            >> $LOG 2>&1

        # Remove temporary files but keep the HMM database files for inspection
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

rule create_plot:
    input:
        contig_file="analysis_other/assembly_info/{asm}/contigs.txt",
        fai_file="analysis_other/assembly_info/{asm}/assembly.fasta.fai",
        feature_bed="analysis_other/features/{asm}/features.bed"
    output:
        plot="analysis_other/plots/{asm}/{asm}_genomic_features.pdf"
    conda:
        "../env/R.yml"
    log:
        "logs/create_plot/{asm}_create_plot.log"
    shell:
        """
        mkdir -p $(dirname {output.plot})

        Rscript workflow/scripts/20_plot_genomic_features.R \
            --contig_file {input.contig_file} \
            --fai_file {input.fai_file} \
            --feature_files {input.feature_bed} \
            --output_prefix analysis_other/plots/{wildcards.asm}/{wildcards.asm} \
            > {log} 2>&1
        """

# Cenmap
rule cenmap:
    input:
        asm = "assembly/output/verkko/{asm}/assembly.fasta",
        mod = lambda wc: find_input_datasets(SimpleNamespace(dataset=wc.asm, type="UL"))["files"][0],
        hq = "assembly/input/{asm}/{asm}.HQ_herro.50x.fastq.gz"
    output:
        done = "analysis_other/cenmap/{asm}/cenmap.done"
    conda:
        "../env/cenmap.yml"
    log:
        "logs/cenmap/{asm}.yml"
    benchmark:
        "runtimes/cenmap/{asm}/cenmap.txt"
    threads:
        24
    shell:
        """
        FA=$(realpath {input.asm})
        HQ=$(realpath {input.hq})
        MOD=$(realpath {input.mod})
        LOG=$(realpath {log})
        WD=$(dirname {output.done})

        pushd $WD >{log}

        cenmap \
            -i $FA \
            -s {wildcards.asm} \
            --hifi $HQ \
            --ont $MOD \
            >$LOG 2>&1
        touch cenmap.done
        """
