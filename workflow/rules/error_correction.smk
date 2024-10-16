#rule preprocess_simplex:
#   input:
#       trim_len =
#       infile = 
#       min_len = 
#   output:
#       ".fastq"
#   conda:#
#
#   shell:
#       """
#       python3 scripts/fq_prep.py -i {input.infile} -o {output} -t {input.trim_len} -l {input.min_len}
#       """
def get_reads_for_sample(wc):
    if wc.sample == "published":
        return{'fastq' : expand("downloads/simplex/{path}.fastq.gz", path = wc.path)}
    

rule dorado_trim:
    input:
        unpack(get_reads_for_sample)
    output:
        fastq = "corrected/{sample}/{path}.trimmed.fastq"
    params: 
        dorado = config['dorado']
    threads:
        60
    benchmark:
        "runtimes/dorado_trim_{sample}_{path}.txt"
    log:
        "logs/dorado_trim_{sample}_{path}.log"
    shell:
        """
        {params.dorado} trim \
            --threads {threads} \
            --emit-fastq \
            {input.fastq} > {output.fastq} \
            2> {log}
        """

rule dorado_correct_mapping:
    input:
        fastq = rules.dorado_trim.output.fastq
    output:
        paf = "corrected/{sample}/{path}.ovl.paf"
    params: 
        dorado = config['dorado'],
        herro_model = config['herro_model']
    threads:
        60
    benchmark:
        "runtimes/dorado_correct_mapping_{sample}_{path}.txt"
    log:
        "logs/dorado_correct_mapping_{sample}_{path}.log"
    shell:
        """
        {params.dorado} correct \
            --to-paf \
            --threads {threads} \
            --model-path {params.herro_model} \
            --device 'cpu' \
            {input.fastq} > {output.paf} \
            2> {log}
        """

rule dorado_correct_inference:
    input:
        fastq = rules.dorado_trim.output.fastq,
        paf = rules.dorado_correct_mapping.output.paf
    output:
        fa = "corrected/{sample}/{path}.corrected.fasta"
    params:
        dorado = config['dorado'],
        herro_model = config['herro_model']
    threads:
        6
    benchmark:
        "runtimes/dorado_correct_inference_{sample}_{path}.txt"
    log:
        "logs/dorado_correct_inference_{sample}_{path}.log"
    resources:
        gpu = 2
    shell:
       """
        {params.dorado} correct \
            --from-paf {input.paf} \
            --threads {threads} \
            --model-path {params.herro_model} \
            --device 'cuda:all' \
            {input.fastq} > {output.fa} \
            2> {log}
       """