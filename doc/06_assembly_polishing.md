# Assembly polishing using Medaka and APK

Assembly polishing with Medaka using a combination of APK (Accuracy Polishing Kit) and ultra-long reads was benchmarked as part of this study. Polishing is integrated into the Snakemake workflow (`workflow/rules/polishing.smk`).

## Medaka command

The polishing step runs `medaka_consensus_joint` with joint APK and ULK models:

```bash
medaka_consensus_joint \
    -i <apk.fastq.gz> -v apk \
    -i <UL.50x.fastq.gz> -v ulk \
    -t 30 \
    -o <output_dir> \
    -d <assembly.fasta> \
    -m r1041_e82_260bps_joint_apk_ulk_v5.0.0 
```

## Runtime

Polishing requires approximately 2 days on 30 threads: roughly 1 day each for the mapping and inference steps.
