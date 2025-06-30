def update_herro_paths(f, dataset):
    f = os.path.basename(f)
    match = re.search(
        r"(\.fastq|\.fastq.gz|\.bam|\.fasta|\.fq|\.fq.gz|\.fa|\.cram)$", f, flags=re.I
    )
    if bool(match):
        ext = match[1]
        fn = re.sub(rf"\{ext}$", "", f)
    else:
        ext = ""
        fn = f
    # print(fn, ext)
    f = os.path.join("data", "corrected", fn, fn + ".corrected.fasta")
    return f


def find_input_datasets(wc, print_debug=True, all_porec=False):
    # Input is a raw data folder, run basecalling depending on type
    # print(datasets[wc.dataset][wc.type])
    files = []
    folders = []
    elements = ()

    from pprint import pprint

    if print_debug:
        print(f"{wc.dataset} | {wc.type}")

    # Handle combined datasets differently
    if str(wc.type).startswith("HQ_combined"):
        (cov_duplex, cov_herro) = str(wc.type).replace("HQ_combined.", "").split("_")
        wc_type = "HQ_combined"
        if print_debug:
            print("--> split combined dataset")
        # Take all duplex reads if called max
        if cov_duplex == "maxx":
            files.extend(
                [f"assembly/input/{wc.dataset}/{wc.dataset}.HQ_duplex.fastq.gz"]
            )
        else:
            files.extend(
                [
                    f"assembly/input/{wc.dataset}/{wc.dataset}.HQ_duplex.{cov_duplex}.fastq.gz"
                ]
            )

        if cov_herro == "maxx":
            files.extend(
                [f"assembly/input/{wc.dataset}/{wc.dataset}.HQ_herro.fastq.gz"]
            )
        else:
            files.extend(
                [
                    f"assembly/input/{wc.dataset}/{wc.dataset}.HQ_herro.{cov_herro}.fastq.gz"
                ]
            )

        elements = {
            "HQ_herro": datasets[wc.dataset]["HQ_herro"]
            if isinstance(datasets[wc.dataset]["HQ_herro"], list)
            else [datasets[wc.dataset]["HQ_herro"]],
            "HQ_duplex": datasets[wc.dataset]["HQ_duplex"]
            if isinstance(datasets[wc.dataset]["HQ_duplex"], list)
            else [datasets[wc.dataset]["HQ_duplex"]],
        }
    else:

        real_type = "POREC" if wc.type == "POREC_all" else wc.type

        elements = {
            wc.type: datasets[wc.dataset].get(real_type, "")
            if isinstance(datasets[wc.dataset].get(real_type, ""), list)
            else [datasets[wc.dataset].get( real_type, "")]
        }

        for wc_type, elements in elements.items():
            # print("Input elements:")
            # pprint(elements)
            for e in elements:
                # Is a folder
                if os.path.isdir(e):
                    if print_debug:
                        print(f"--> transform input folder: {e}")

                    if wc_type == "UL":
                        files.append(
                            f"data/basecalled/sup/{os.path.basename(e)}/{os.path.basename(e)}.sup.unmapped.bam"
                        )
                    if wc_type == "HQ_duplex":
                        files.append(
                            f"data/basecalled/Duplex/{os.path.basename(e)}/{os.path.basename(e)}.duplexonly.unmapped.bam"
                        )
                    if wc_type == "HQ_herro":
                        files.append(
                            f"data/corrected/{os.path.basename(e)}/{os.path.basename(e)}.corrected.fasta"
                        )
                    if wc_type == "HQ_paternal":
                        files.append(
                            f"data/corrected/{os.path.basename(e)}/{os.path.basename(e)}.corrected.fasta"
                        )
                    if wc_type == "HQ_maternal":
                        files.append(
                            f"data/corrected/{os.path.basename(e)}/{os.path.basename(e)}.corrected.fasta"
                        )
                    if wc_type == "POREC":
                        files.append(
                            f"data/basecalled/sup/{os.path.basename(e)}/{os.path.basename(e)}.sup.unmapped.bam"
                        )
                        if not all_porec: break # Only use first POREC for assembly
                    if wc_type == "POREC_all":
                        files.append(
                            f"data/basecalled/sup/{os.path.basename(e)}/{os.path.basename(e)}.sup.unmapped.bam"
                        )
                    if wc_type == "APK":
                        files.append(
                            f"data/basecalled/apk/{os.path.basename(e)}/{os.path.basename(e)}.apk.unmapped.bam"
                        )
                # Is a file
                else:
                    match wc_type:
                        case "HQ_herro":
                            f = update_herro_paths(e, wc.dataset)
                            if print_debug:
                                print(f"--> transform herro path: {e}")
                            files.append(f)
                        case "HQ_paternal":
                            f = update_herro_paths(e, wc.dataset)
                            if print_debug:
                                print(f"--> transform paternal reads path: {e}")
                            files.append(f)
                        case "HQ_maternal":
                            f = update_herro_paths(e, wc.dataset)
                            if print_debug:
                                print(f"--> transform maternal reads path: {e}")
                            files.append(f)
                        case "HQ_duplex":
                            files.append(e)
                        case "UL":
                            files.append(e)
                        case "POREC":
                            files.append(e)
                            if not all_porec: break # Only use first POREC for assembly
                        case "POREC_all":
                            files.append(e)
                        case "APK":
                            files.append(e)
                        case _:
                            if print_debug:
                                print(f"Unrecognized dataset type {wc.type} for {e}")
                            files.append(e)

    if not files:
        if print_debug:
            print("Files empty: files generated by other rule")
    else:
        if print_debug:
            print("Updated files:")
        if print_debug:
            pprint(files)
    if print_debug:
        print("-" * 100)

    return {
        "files": files,
        "folders": folders,
    }


def input_isbam(wc):
    files = find_input_datasets(wc, print_debug=False)["files"]
    if all(f.endswith(".bam") for f in files):
        return True
    elif all(not f.endswith(".bam") for f in files):
        return False
    else:
        raise ValueError(
            "Inconsistent file types: some files have .bam ending and some don't. Cannot merge"
        )


rule merge_copy_rename_fastq:
    input:
        unpack(find_input_datasets),
    output:
        "assembly/input/{dataset}/{dataset}.{type}.fastq.gz",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/merge_copy_rename_fastq/{dataset}.{type}.log",
    params:
        samtools=lambda wc: "samtools" if input_isbam(wc) else "",
    threads: 12
    shell:
        """
        samtools fastq <({params.samtools} cat {input.files}) 2>{log}\
        | pigz -p {threads} -c >{output} 2>>{log}
        """

rule map_unaligned_bam:
    input:
        bam="data/{path}.bam",
        ref=config["ref"],
    output:
        bam="data/mapped/{path}.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_bam/{path}.log",
    threads: 60
    params:
        preset=lambda wc: "lr:hq" if "duplex" in wc.path else "map-ont",
    shell:
        """
        samtools fastq {input.bam} \
        | minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.ref} - 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

        samtools index {output.bam}
        """


rule map_fq:
    input:
        fq="assembly/input/{file}.fastq.gz",
        ref=config["ref"],
    output:
        bam="data/mapped/{file}.bam",
    conda:
        "../env/minimap2.yml"
    log:
        "logs/map_fq/{file}.log",
    threads: 60
    params:
        preset=lambda wc: "lr:hq" if "duplex" in wc.file else "map-ont",
    shell:
        """
        minimap2 --MD -ax {params.preset} --eqx \
            -t {threads} \
            {input.ref} {input.fq} 2>{log} \
        | samtools sort -m 4G -@ 4 -o {output.bam} -O BAM - >>{log} 2>&1

        samtools index {output.bam}
        """


rule sample_to_target_cov:
    input:
        fq="assembly/input/{file}.fastq.gz",
    output:
        fq=r"assembly/input/{file}.{cov,\d+x}.fastq.gz",
    conda:
        "../env/filtlong.yml"
    log:
        "logs/filtlong/{file}_{cov}.log",
    params:
        min_length=lambda wc: config["min_length"]["UL"] if (".UL" in str(wc.file)) else config["min_length"]["default"],
        min_mean_q=config["min_mean_q"],
        target_base=lambda wc: str(
            int(wc.cov.replace("x", "")) * config["genome_length"]
        ),
    shell:
        """
        filtlong \
            --min_length {params.min_length} \
            --target_bases {params.target_base} \
            --length_weight 10 \
            --min_mean_q {params.min_mean_q} \
            {input} 2>{log} \
        | gzip -c > {output} 2>>{log}
        """


rule extract_location_data:
    input:
        bam="data/mapped/{file}.bam",
    output:
        fq="assembly/input/{file}.{roi,chr.*}.fastq.gz",
    log:
        "logs/extract_location_data/{file}_{roi}.log",
    shell:
        """
        samtools view -h {input.bam} {wildcards.roi} 2>>{log} \
        | samtools fastq 2>>{log} \
        | gzip -c > {output.fq} 2>>{log}
        """


ruleorder: extract_location_data > merge_copy_rename_fastq
