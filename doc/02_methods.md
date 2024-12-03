# Methods

Analysis methods. Collection of informations in instructions for tools used in the pipeline. Comparison and evualuation of different analysis approaches. Small code snippets.


### Basecalling

Basecalling is not yet integrated into the pipeline. Should be done when all samples are assembled in routine.
*Include MOD basecalling for UL reads in next basecalling*

#### #04339
```
bin/dorado-0.8.3-linux-x64/bin/dorado basecaller \
    data/ref/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    data/raw/set1/24070/24070LRa003_04399/20241111_1425_2A_PBA32113_7de81dc2 \
    --recursive \
    --trim all \
    --output-dir data/basecalled/24070LRa003_04399

bin/dorado-0.8.3-linux-x64/bin/dorado summary \
    data/basecalled/24070LRa003_04399/*.bam \
    > data/basecalled/24070LRa003_04399/sequencing_summary.txt
```

#### 04400
```
bin/dorado-0.8.3-linux-x64/bin/dorado basecaller \
    data/ref/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    /mnt/storage3/raw_data/VULCAN/24070/24070LRa010_04400/20241111_1433_2B_PBA36199_3fb5db31 \
    --recursive \
    --trim all \
    --output-dir data/basecalled/24070LRa010_04400 
```
