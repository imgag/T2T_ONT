# Methods

Analysis methods. Collection of informations in instructions for tools used in the pipeline. Comparison and evualuation of different analysis approaches. Small code snippets.

## Use Herro corrected reads for assembly

[Reply](https://github.com/marbl/verkko/issues/322) from S.Koren about best practive for Herro input to T2T 

> You can see some info on the README as well as #321. We usually use longest 50x or whatever can be corrected in a reasonable time with Herro. The suggestions pertaining to minimum k-mer coverage in #321 can be used to address the increasing complexity in the graph with high coverage. There are uncorrected regions with systematic errors that start confirming each other and need to be filtered. You can also try the hifiasm correction as an alternative to herro.
> You should not input corrected reads as UL data. The correction can trim and lose parts of the sequence so you want to allow verkko to fill those gaps with the original UL data. Read length is fine, longer reads will build a more resolved initial graph.

## Strategies described for assembly polishing:

https://github.com/arangrhie/T2T-Polish/tree/master


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
