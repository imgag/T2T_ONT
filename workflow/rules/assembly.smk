rule verkko:
    input:
        ul = "data/{dataset}/{dataset}.UL.fastq.gz",
        hq = "data/{dataset}/{dataset}.HQ.fastq.gz",
        porec = "data/{dataset}/{dataset}.POREC.fastq.gz"
    output:
        gfa = "assembly/{dataset}/assembly.homopolymer-compressed.noseq.gfa",
        fa = "assembly/{dataset}/assembly.fasta",
        hp1 = "assembly/{dataset}/assembly.haplotype1.fasta",
        hp2 = "assembly/{dataset}/assembly.haplotype2.fasta"
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
