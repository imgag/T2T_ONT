# Data QC

Documentation of sequencing quality control for the study dataset.

## Sequencing overview

Flowcells used for the publication dataset:

| Type | Count |
|:-----|------:|
| APK | 3 |
| LSK | 4 |
| PoreC | 49 |
| UltraLong | 63 |
| **total** | 119 |

## Comparison with published UL data (HG002)

We compared our first UL sequencing runs against the published UL dataset from Koren et al. (2024).

### Characteristics of the published dataset

Sequencing was performed in 2022 on R9.4.1 flowcells across 8 separate flowcells. The available download appears to be pre-filtered for read length.

### Comparison before filtering

| filename     | yield        | N50   | max readlen | yield>100kb | yield>200kb | reads > 1Mb |
| :----------- | :----------- | :---- | :---------- | :---------- | :---------- | :---------- |
| published.UL | 496106934638 | 81748 | 1.9 Mb      | 190 Gb      | 31 Gb       | 45          |
| run_04399.UL | 113914796731 | 67652 | 1.7 Mb      | 38 Gb       | 10 Gb       | 29          |
| run_04400.UL | 107609849976 | 77925 | 2.0 Mb      | 40 Gb       | 9 Gb        | 54          |

Our per-flowcell yield exceeds the average of the published dataset (62.5 Gb/flowcell after filtering). Mean read quality is substantially better with the newer R10.4.1 pores.

### After filtering to 50x coverage

| filename         | yield | N50    | max readlen | yield>100kb | yield>200kb | reads > 1Mb |
| :--------------- | :---- | :----- | :---------- | :---------- | :---------- | :---------- |
| TUE_01.UL.50x    | 160GB | 97739  | 2037kb      | 78GB        | 18GB        | 50          |
| published.UL.50x | 160GB | 129072 | 1511kb      | 121GB       | 22GB        | 39          |

At matched coverage (50x), the published R9.4.1 dataset has higher N50 and greater ultra-long yield (>100 kb), reflecting the maturity of that protocol. Our R10.4.1 data achieves better per-read quality.

## Coverage requirements

We need at least 50x long read coverage (~160 Gb) of ultra-long reads. Benchmarking showed that 50x coverage yields similar assembly performance to the 70x commonly recommended in the literature.

## Downsampling strategy

To benchmark required read depths for T2T assembly, the published HG002 dataset was downsampled using `filtlong`, preferring high-quality and longer reads. Mapping coverages were calculated to confirm that target depths were achieved.

![Mapping coverages for benchmark assembly input](img/subsampling_coverages.cov_boxplot.png)
