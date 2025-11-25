rule all_extended_qc:
    input:        
        # CenMap
        expand("analysis_other/cenmap/{asm}/cenmap.done", asm = finished_samples)

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