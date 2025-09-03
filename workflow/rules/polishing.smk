# This path is for polishing the final assembly using APK and UL data
rule medaka_polishing:
    input:
        ul = lambda wc: get_assembly_input(wc).get('ul'),
        apk = lambda wc: get_assembly_input(wc).get('apk'),
        asm = "assembly/output/verkko/{asm}/assembly.fasta"
    output:
        "assembly/output/verkko/{asm}/assembly_polished.fasta"
    log:
        "logs/medaka_polishing/{asm}.log"
    benchmark:
        "runtimes/{asm}.medaka_polishing.txt"
    threads: 60
    conda:
        "../env/medaka.yml"
    shell:
        """
        medaka_consensus_joint \
            -i {input.apk} -v apk \
            -i {input.ul} -v ulk \
            -d {input.asm} \
            -t {threads} \
            -m r1041_e82_260bps_joint_apk_ulk_v5.0.0 \
            -o $(dirname {output})/polishing \
            >{log} 2>&1
        mv -v $(dirname {output})/polishing/consensus.fasta {output} >> {log} 2>&1
        rm -rfv $(dirname {output})/polishing >>{log} 2>&1
        """

# possible todo add polishing with dorado and UL?