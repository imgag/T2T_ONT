def get_assembly_input(wc):
    s_d = asm[wc.asm]['dataset']
    s_hq = asm[wc.asm]['HQ_method']
    s_hq = "HQ_" + s_hq
    s_cov_ul = asm[wc.asm].get('cov_UL', '')
    if s_cov_ul != '':
        s_cov_ul = str(s_cov_ul)+'x.'
    s_cov_hq = asm[wc.asm].get('cov_HQ', '')
    if s_cov_hq != '':
        s_cov_hq = str(s_cov_hq)+'x.'
    s_re = asm[wc.asm].get('region', '')
    if s_re != '':
        s_re = s_re+'.'
#    print(s_d, s_ul, s_hq, s_re)
    files = {
        'ul' : f"assembly/input/{s_d}/{s_d}.UL.{s_cov_ul}{s_re}fastq.gz",
        'hq' : f"assembly/input/{s_d}/{s_d}.{s_hq}.{s_cov_hq}{s_re}fastq.gz",
        'porec' : f"assembly/input/{s_d}/{s_d}.POREC.fastq.gz"
    }
    return(files)

rule verkko:
    input:
        unpack(get_assembly_input)
    output:
        gfa = "assembly/output/{asm}/assembly.homopolymer-compressed.noseq.gfa",
        fa = "assembly/output/{asm}/assembly.fasta",
        hp1 = "assembly/output/{asm}/assembly.haplotype1.fasta",
        hp2 = "assembly/output/{asm}/assembly.haplotype2.fasta"
    conda:
        "../env/verkko.yml"
    log:
        "logs/verkko_{asm}.log"
    benchmark:
        "runtimes/{asm}.verkko.txt"
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
