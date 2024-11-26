def get_assembly_input(wc):
    s_d = asm[wc.asm]["dataset"]
    s_hq = asm[wc.asm]["HQ_method"]
    s_hq = "HQ_" + s_hq
    s_cov_ul = asm[wc.asm].get("cov_UL", "")
    if s_cov_ul != "":
        s_cov_ul = str(s_cov_ul) + "x."
    s_cov_hq = asm[wc.asm].get("cov_HQ", "")
    if s_cov_hq != "":
        s_cov_hq = str(s_cov_hq) + "x."
    s_re = asm[wc.asm].get("region", "")
    if s_re != "":
        s_re = s_re + "."
    files = {
        "ul": f"assembly/input/{s_d}/{s_d}.UL.{s_cov_ul}{s_re}fastq.gz",
        "hq": f"assembly/input/{s_d}/{s_d}.{s_hq}.{s_cov_hq}{s_re}fastq.gz",
        "porec": f"assembly/input/{s_d}/{s_d}.POREC.fastq.gz",
    }
    return files


rule verkko:
    input:
        unpack(get_assembly_input),
    output:
        gfa=temp("/tmp/verrko_{asm}/assembly.homopolymer-compressed.noseq.gfa"),
        fa=temp("/tmp/verrko_{asm}/assembly.fasta"),
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/verkko_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko.txt"
    threads: 320
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
    shell:
        """
        verkko -d $(dirname {output.fa}) \
            --hifi {input.hq} \
            --nano {input.ul} \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """


# Second rule for scaffolding is needed because not all required outputs are well defined in verrko. This workaround helps recreate all required output files if they are deleted.
rule verkko_scaffold:
    input:
        unpack(get_assembly_input),
    output:
        hp1=temp("/tmp/verrko_{asm}/assembly.haplotype1.fasta"),
        hp2=temp("/tmp/verrko_{asm}/assembly.haplotype2.fasta"),
        colors=temp("/tmp/verrko_{asm}/assembly.colors.csv"),
    conda:
        "../env/verkko.yml"
    group:
        "verrko"
    log:
        "logs/verkko_scaffold_{asm}.log",
    benchmark:
        "runtimes/{asm}.verkko_scaffold.txt"
    threads: 120
    params:
        dryrun="--dryrun" if config["verkko_dryrun"] else "",
    shell:
        """
        verkko -d $(dirname {output.hp1}) \
            --hifi {input.hq} \
            --nano {input.ul} \
            --porec {input.porec} \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """


rule verkko_copy_results:
    input:
        gfa="/tmp/verrko_{asm}/assembly.homopolymer-compressed.noseq.gfa",
        fa="/tmp/verrko_{asm}/assembly.fasta",
        hp1="/tmp/verrko_{asm}/assembly.haplotype1.fasta",
        hp2="/tmp/verrko_{asm}/assembly.haplotype2.fasta",
        colors="/tmp/verrko_{asm}/assembly.colors.csv",
    output:
        gfa="assembly/output/{asm}/assembly.homopolymer-compressed.noseq.gfa",
        fa="assembly/output/{asm}/assembly.fasta",
        hp1="assembly/output/{asm}/assembly.haplotype1.fasta",
        hp2="assembly/output/{asm}/assembly.haplotype2.fasta",
        colors="assembly/output/{asm}/assembly.colors.csv",
    log:
        "logs/verkko_copy_results_{asm}.log",
    group:
        "verrko"
    shell:
        """
        cp -v /tmp/verrko_{wildcards.asm}/assembly.* assembly/output/{wildcards.asm}/  >{log} 2>{log}
        rm -rf /tmp/verrko_{wildcards.asm}
        """
