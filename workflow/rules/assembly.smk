rule verkko:
    input:
        ul = "assembly/input/{dataset}/{dataset}.UL.fastq.gz",
        hq = "assenmbly/input/{dataset}/{dataset}.HQ.fastq.gz",
        porec = "assembly/input/{dataset}/{dataset}.POREC.fastq.gz"
    output:
        gfa = "assembly/output/{dataset}/assembly.homopolymer-compressed.noseq.gfa",
        fa = "assembly/output/{dataset}/assembly.fasta",
        hp1 = "assembly/output/{dataset}/assembly.haplotype1.fasta",
        hp2 = "assembly/output/{dataset}/assembly.haplotype2.fasta"
    conda:
        "../env/verkko.yml"
    log:
        "logs/verkko_{dataset}.log"
    benchmark:
        "runtimes/{dataset}.verkko.txt"
    threads:
        120
    params:
        dryrun = "--dryrun" if config['verkko_dryrun'] else ""
    shell:
        """
        verkko -d $(dirname {output.fa}) \
            --hifi {input.hq} \
            --nano {input.ul} \
            --porec {input.porec} \
            --snakeopts "--cores {threads} {params.dryrun}" \
            >{log} 2>{log}
        """
