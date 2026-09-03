rule get_contig_stats:
    input:
        unpack(get_assembly_output),
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/contig_stats.{hp}.tsv",
    log:
        "logs/get_contig_stats/{isphased}_{tool}_{asm}_{hp}.tsv",
    conda:
        "../env/samtools.yml"
    threads: 1
    shell:
        """
        samtools faidx {input.assembly}
        awk 'BEGIN {{ count_total=0; count_gt_10Mb=0; }} {{ current_length=$2+0; count_total++; if (current_length > 10000000) count_gt_10Mb++; }} END {{ printf "n_contigs\\t%d\\n", count_total; printf "n_contigs_over_10mb\\t%d\\n", count_gt_10Mb; }}' \
        {input.assembly}.fai > {output} 2>{log}
        """


rule subsample_ref_genome:
    input:
        fa="data/ref/{ref}.fasta",
    output:
        fa="data/ref/{ref}.{roi,chr.*}.fasta",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/subsample_ref_genome/{ref}_{roi}.log",
    threads: 1
    shell:
        """
        samtools faidx {input.fa} 2>{log}
        samtools faidx {input.fa} {wildcards.roi} > {output} 2>>{log}
        samtools faidx {output.fa} 2>>{log}
        """


rule create_colors:
    input:
        paf="assembly/qc/{isphased}_{tool}/{asm}/both.mapped_T2T.paf",
        colors_phasing=lambda wc: (
            f"assembly/output/verkko/{wc.asm}/assembly.colors.csv"
            if wc.isphased == "phased" and wc.tool == "verkko"
            else []
        ),
    output:
        csv="assembly/qc/{isphased}_{tool}/{asm}/colors.tsv",
    log:
        "logs/create_colors/{isphased}_{tool}_{asm}.log",
    params:
        colours_phasing=lambda wc: (
            f"-c assembly/output/verkko/{wc.asm}/assembly.colors.csv"
            if wc.isphased == "phased" and wc.tool == "verkko"
            else ""
        ),
    shell:
        """
        python workflow/scripts/11_extract_colors.py \
            -i {input.paf} {params.colours_phasing} \
            -o {output.csv} \
            >{log} 2>&1
        """


rule process_graph:
    input:
        gfa=get_assembly_graph_output,
        scfmap="assembly/output/verkko/{asm}/assembly.scfmap",
        color="assembly/qc/{isphased}_{tool}/{asm}/colors.tsv",
    output:
        gfa="assembly/qc/{isphased}_{tool}/{asm}/assembly_graph.gfa",
    log:
        "logs/process_graph/{isphased}_{tool}_{asm}.log",
    shell:
        """
        python workflow/scripts/12_process_gfa.py \
            --gfa {input.gfa} \
            --scfmap {input.scfmap} \
            --colors {input.color} \
            --output {output.gfa} \
            > {log} 2>&1
        """


rule bandage:
    input:
        gfa="assembly/qc/{isphased}_{tool}/{asm}/assembly_graph.gfa",
        color=get_assembly_colors,
    output:
        svg="assembly/qc/{isphased}_{tool}/{asm}/bandage_graph.svg",
        png="assembly/qc/{isphased}_{tool}/{asm}/bandage_graph.png",
    conda:
        "../env/bandage.yml"
    log:
        "logs/bandage/{isphased}_{tool}_{asm}.log",
    threads: 1
    shell:
        """
        BandageNG image {input.gfa} {output.svg} --color {input.color} > {log} 2>&1
        BandageNG image {input.gfa} {output.png} --color {input.color} > {log} 2>&1
        """


rule qc_paftools_stat:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        ref=get_ref_genome,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_stat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_stat/{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js stat\
            {input.paf} \
            > {output} 2>{log}
        """


rule qc_paftools_asmstat:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        ref=get_ref_genome,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_asmstat.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmstat/{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmstat\
            {input.ref}.fai {input.paf} \
            > {output} 2>{log}
        """


rule qc_paftools_misjoin:
    input:
        paf=rules.map_asm_to_ref.output.paf,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_misjoin.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_misjoin/{isphased}_{tool}_{asm}_{hp}.log",
    params:
        centromere_bed=config["ref_centromer_bed"],
    shell:
        """
        paftools.js misjoin\
            -c {params.centromere_bed} \
            -e \
            {input.paf} \
            > {output} 2>{log}
        """


def get_ref_cdna_paf(wc):
    import re

    ref = config["ref"].replace(".fasta", ".cdna.paf")
    # print(ref)
    match = re.search(r"chr\d+", str(wc.asm))
    # print(match)
    if match:
        # print(match)
        prefix = config["ref"].replace("fasta", "")
        ref = f"{prefix}{match.group(0)}.cdna.paf"
    return ref


rule qc_paftools_asmgene:
    input:
        paf_asm=rules.map_cdna_to_asm.output.paf,
        paf_ref=get_ref_cdna_paf,
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/qc_paftools_asmgene.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/paftools_asmgene/{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        paftools.js asmgene \
            {input.paf_ref} {input.paf_asm} \
            > {output} 2>{log}
        """


rule scaffold_lengths:
    input:
        unpack(get_assembly_output),
        ref=get_ref_genome,
    output:
        txt="assembly/qc/phased_{tool}/{asm}/scaffold_lengths.{hp}.txt",
    conda:
        "../env/minimap2.yml"
    threads: 1
    log:
        "logs/scaffold_lengths/{tool}_{asm}_{hp}.txt",
    shell:
        """
        samtools faidx {input.assembly}
        cut -f 1,2 {input.assembly}.fai > {output}
        cut -f 1,2 {input.ref}.fai >> {output}
        rm {input.assembly}.fai
        """

def get_done_file(wc):
    if wc.isphased == "phased" and wc.tool == "hprc":
        return []
    elif wc.isphased == "phased" and wc.tool == "verkko":
        return f"assembly/output/verkko/{wc.asm}/create_scaffold.done"
    elif asm[wc.asm].get("trio_phasing", False):
        return f"assembly/output/verkko/{wc.asm}/create_scaffold.done"
    else:
        return []

rule gap_stats:
    input:
        unpack(get_assembly_output),
        done=get_done_file,
    output:
        bed="assembly/qc/{isphased}_{tool}/{asm}/gap_stats.{hp}.n_regions.bed",
        stats="assembly/qc/{isphased}_{tool}/{asm}/gap_stats.{hp}.n_stats.tsv",
    conda:
        "../env/hapdiff.yml"
    log:
        "logs/gap_stats/{isphased}_{tool}_{asm}_{hp}.log",
    shell:
        """
        workflow/scripts/17_assembly_qc_gaps.py \
            --input {input.assembly} \
            --output $(dirname {output.bed})/gap_stats.{wildcards.hp} \
            > {log} 2>&1
        """


rule qc_meryl:
    input:
        ref_q100=config["ref_hg002_q100"],
    output:
        meryl=directory("data/ref/hg002_q100_meryl/hg002_q100_k_{k_val}.meryl"),
    params:
        k=config["K-mer"],
    conda:
        "../env/merqury.yml"
    log:
        "logs/meryl_count/ref_q100_{k_val}.log",
    shell:
        """
        meryl count k={params.k} {input.ref_q100} output {output.meryl} > {log} 2>&1
        """

rule qc_meryl_shortread:
    input:
        fq = lambda wc: config["kmer_refs"].get(wc.asm, ""),
    output:
        meryl=directory("analysis_other/merqury_shortread/{asm}/{asm}_{k_val}.meryl"),
    params:
        k=config["K-mer"],
    conda:
        "../env/merqury.yml"
    resources:
        mem_gb=500,
    log:
        "logs/meryl_count/{asm}_{k_val}.log",
    shell:
        """
        meryl count k={params.k} {input.fq} output {output.meryl} > {log} 2>&1
        """

# Select correct Meryl ref for sample. Two options:
# 1) HG002 for published set, fallback 2) Illumina shortread data from config file for specified samples
def get_meryl_ref(wc):
    if "published" in str(wc.asm):
        return f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl'
    elif str(wc.asm) in config["kmer_refs"]:
        ref = config["kmer_refs"].get(wc.asm, "")
        if ref.endswith(".meryl"):
            return ref
        else:
            return(f'analysis_other/merqury_shortread/{wc.asm}/{wc.asm}_k_{config["K-mer"]}.meryl')

    elif "TUE_02" in str(wc.asm):  # Keep for old versions of T2T00 samplename
        return "analysis_other/merqury_shortread/DX203429_02.meryl"
    else:
        print(f"WARNING: Unknown Meryl Ref for {wc.asm}, using HG002")
        return f'data/ref/hg002_q100_meryl/hg002_q100_k_{config["K-mer"]}.meryl'


rule qc_merqury_haplotypes:
    input:
        meryl=get_meryl_ref,
        pat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype1"})["assembly"],
        mat_fa=lambda wc: get_assembly_output({**wc, "hp": "haplotype2"})["assembly"],
    output:
        out="assembly/qc/phased_{tool}/{asm}/merqury.qv",
        hap_pat_meryl="assembly/qc/phased_{tool}/{asm}/merqury.haplotype1.qv",
        hap_mat_meryl="assembly/qc/phased_{tool}/{asm}/merqury.haplotype2.qv",
    conda:
        "../env/merqury.yml"
    log:
        "logs/merqury/{tool}_{asm}.log",
    threads: 40
    shell:
        """
        INPUT_MERYL=$(realpath {input.meryl})
        INPUT_PAT_FA=$(realpath {input.pat_fa})
        INPUT_MAT_FA=$(realpath {input.mat_fa})
        OUTPUT_PREFIX=$(dirname $(realpath {output.out}))/merqury
        LOG_FILE=$(realpath {log})
        pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
        cp $INPUT_PAT_FA haplotype1.fa
        cp $INPUT_MAT_FA haplotype2.fa
        export PATH=$PATH:"$CONDA_PREFIX"/share/merqury/eval
        qv.sh \
            $INPUT_MERYL \
            haplotype1.fa \
            haplotype2.fa \
            $OUTPUT_PREFIX \
         >> $LOG_FILE 2>&1
        
        rm -r *.meryl >> $LOG_FILE 2>&1
        popd >> $LOG_FILE 2>&1
        """


def get_merqury_input(wc):
    if wc["tool"] == "verkko":
        if wc["isphased"] == "unphased":
            return {'fa': f"assembly/output/verkko_unphased/{wc['asm']}/assembly{wc['polished']}.fasta"}
        elif wc["isphased"] == "phased":
            # TODO this does not work for polished trios, not needed right now
            return {
                "fa": f"assembly/output/verkko/{wc['asm']}/assembly{wc['polished']}.fasta",
                "done": f"assembly/output/verkko/{wc['asm']}/create_scaffold.done"
            }
    elif wc["tool"] == "hifiasm":
        return {"fa": f"assembly/output/hifiasm/{wc['asm']}/assembly{wc['polished']}.fasta"}
    elif wc["tool"] == "repair":
        return {"fa": f"analysis_other/assembly_repairer/{wc['asm']}/repaired_assembly.final.fa"}
    else :
        raise ValueError(f"Merqury input: Invalid tool{wc['tool']}")

rule qc_merqury_both:
    input:
        unpack(get_merqury_input),
        meryl=get_meryl_ref,
    output:
        out="assembly/qc/{isphased}_{tool}/{asm}/merqury{polished}.qv",
    wildcard_constraints:
        polished=".polished*",
    conda:
        "../env/merqury.yml"
    log:
        "logs/merqury/{isphased}_{tool}_{asm}{polished}.log",
    threads: 40
    shell:
        """
        INPUT_MERYL=$(realpath {input.meryl})
        INPUT_FA=$(realpath {input.fa})
        OUTPUT_PREFIX=$(dirname $(realpath {output.out}))/merqury{wildcards.polished}
        LOG_FILE=$(realpath {log})
        pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
        export PATH=$PATH:"$CONDA_PREFIX"/share/merqury/eval
        qv.sh \
            $INPUT_MERYL \
            $INPUT_FA \
            $OUTPUT_PREFIX \
         >> $LOG_FILE 2>&1
        rm -r *.meryl >> $LOG_FILE 2>&1
        popd >> $LOG_FILE 2>&1
        """


ruleorder: qc_merqury_haplotypes > qc_merqury_both


# temporary meryl files are created during qv.sh process so i put the rm
# export PATH=$PATH:"$CONDA_FREFIX"/share/merqury/eval
# because the qv script

rule findt2t_extract_telomere_regions:
    input:
        unpack(get_assembly_output),
    output:
        fa=temp("assembly/qc/{isphased}_{tool}/{asm}/telomere_regions.{hp}.fa"),
    log:
        "logs/findt2t/extract_telomere_regions/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 6
    params:
        seqkit = config.get("seqkit_path", "seqkit"),
    shell:
        """
        # Extract first 1000bp (marked as _Start)
        {params.seqkit} subseq -w 0 -j {threads} -r 1:1000 {input.assembly} | \
            awk '/^>/ {{print $0 "_Start"; next}} {{print}}' > {output.fa} 2>{log}
        
        # Extract last 1000bp (marked as _End)
        {params.seqkit} subseq -w 0 -j {threads} -r -1000:-1 {input.assembly} | \
            awk '/^>/ {{print $0 "_End"; next}} {{print}}' >> {output.fa} 2>>{log}
        """

rule findt2t_extract_seqinfo:
    input:
        assembly=lambda wc: get_assembly_output(wc)["assembly"],
    output:
        seqinfo="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
    log:
        "logs/findt2t/extract_seqinfo/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 6
    params:
        seqkit = config.get("seqkit_path", "seqkit"),
    shell:
        """
        # Extract sequence information
        {params.seqkit} fx2tab \
            -j {threads} \
            -C N -l -n \
            -o {output.seqinfo} \
            {input.assembly} >{log} 2>&1
        """

rule findt2t_detect_telomere_motifs:
    input:
        fa="assembly/qc/{isphased}_{tool}/{asm}/telomere_regions.{hp}.fa",
    output:
        csv="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif.csv",
    log:
        "logs/detect_telomere_motifs/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 6
    params:
        tidk = config.get("tidk_path", "tidk"),
        motif=lambda wc: config.get("telomere_motif", "TTAGGG"),
        output_prefix=lambda wc: f"telomere_windows.{wc.hp}",
    shell:
        """  
        {params.tidk} search -d $(dirname {output.csv}) \
            --fasta {input.fa} \
            --output {params.output_prefix} \
            --string {params.motif} \
            --window 1000 >>{log} 2>&1
        
        mv -v $(dirname {output.csv})/{params.output_prefix}_telomeric_repeat_windows.csv {output.csv} 2>>{log} || touch {output.csv}
        """

rule findt2t_identify_motif_contigs:
    input:
        csv="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif.csv",
        seqinfo="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
    output:
        motif="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif_T2T.txt",
    log:
        "logs/identify_motif_contigs/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 1
    shell:
        """
        # Extract T2T contigs with both telomere ends
        tail -n +2 {input.csv} | \
            awk -F "," '$3 + $4 >=15{{print $1}}' | \
            while IFS= read -r line; do
                contig=${{line%_*}}
                awk -F '\\t' -v id="$contig" '$1==id{{print; exit}}' {input.seqinfo}
            done | \
            cut -f1,2 | sort | uniq -c | sort -rgk1 | \
            tr -s "[:blank:]" "\\t" | sed 's/^\\t//' | \
            {{ grep "^2" || true; }} | cut -f2- > {output.motif} 2>>{log}
        
        # Ensure output file exists even if empty
        touch {output.motif}
        """

rule findt2t_calculate_alignment_coverage:
    input:
        paf=rules.map_asm_to_ref.output.paf,
        seqinfo="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
        motif_csv="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif.csv",
    output:
        alignment="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.alignment.txt",
    log:
        "logs/calculate_alignment_coverage/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 1
    params:
        cov_cal = config.get("covcal_path", "cov_cal"),
    shell:
        """
        # Find continuous alignments with over 0.95 coverage (exits with code 1)
        {params.cov_cal} -T {input.paf} >{log} 2>{output.alignment} || touch {output.alignment}
        """

rule findt2t_filter_t2t_contigs:
    input:
        seqinfo="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.seqinfo.txt",
        alignment_csv="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}.alignment.txt",
        motif_csv="assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_motif.csv",
    output:
        "assembly/qc/{isphased}_{tool}/{asm}/T2T_contigs.{hp}_alignment_T2T.txt",
    log:
        "logs/filter_t2t_contigs/{isphased}_{tool}_{asm}_{hp}.log",
    threads: 1
    shell:
        """
        # Extract T2T contigs with both telomere ends and good alignment
        tail -n +2 {input.motif_csv} | \
            awk -F "," '$3 + $4 >= 15{{print $1}}' | \
            while IFS= read -r line; do
                contig=${{line%_*}}
                seq_record=$(awk -F '\\t' -v id="$contig" '$1==id{{print; exit}}' {input.seqinfo})
                aln_record=$(awk -F '\\t' -v id="$contig" '$1==id{{print; exit}}' {input.alignment_csv})
                if [[ -n "$seq_record" && -n "$aln_record" ]]; then
                    printf '%s\\t%s\\n' "$seq_record" "$aln_record"
                fi
            done | \
            cut -f 1-3,9,10 | awk 'NF>4{{print}}' | \
            sort | uniq -c | tr -s "[:blank:]" "\\t" | \
            sed 's/^\\t//' | {{ grep "^2" || true; }} | cut -f2- > {output} 2>>{log}
        
        # Ensure output file exists even if empty
        touch {output}
        """



rule find_T2T_scaffolds:
    input:
        motif_file="assembly/qc/{isphased}_{tool}/{sample}/T2T_contigs.{hp}_motif_T2T.txt",
        alignment_file="assembly/qc/{isphased}_{tool}/{sample}/T2T_contigs.{hp}_alignment_T2T.txt",
        gap_file="assembly/qc/{isphased}_{tool}/{sample}/gap_stats.{hp}.n_regions.bed"
    output:
        "assembly/qc/{isphased}_{tool}/{sample}/T2T_scaffolds.{hp}.txt"
    log:
        "logs/find_T2T_scaffolds/{isphased}_{tool}_{sample}_{hp}.log",
    shell:
        """
        python3 workflow/scripts/41_collect_T2T_contigs.py \
            --samples {wildcards.sample} \
            --motif-files {input.motif_file} \
            --alignment-files {input.alignment_file} \
            --gap-files {input.gap_file} \
            --output {output} \
            > {log} 2>&1
        """