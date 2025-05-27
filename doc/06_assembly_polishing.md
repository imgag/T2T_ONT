# Assembly polishing using Medaka and APK

Test command to run medaka:

```bash
medaka_consensus_joint \
    -i /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/basecalled/TUE02.apk.fastq.gz -v apk \
    -i /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/input/TUE_02/TUE_02.UL.50x.fastq.gz -v ulk \
    -t 30 \
    -o medaka_test \
    -d /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/output/verkko/TUE_02/assembly.fasta \
    -m r1041_e82_260bps_joint_apk_ulk_v5.0.0 
```

## Runtimes

Started: Wednesday 12.03 ~12:00 (30 Threads)
Ended: Friday 14.03 ~17:30

~ 1 day for each Mapping and Inference.
