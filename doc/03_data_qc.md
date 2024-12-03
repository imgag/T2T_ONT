# Data QC
Comparison of our seq data to public datasets. Validation of correct data preprocessing. Sequencing recommendations.

## Ultra Long ONT Reads

In this document we compare our first UL libraries with published UL dataset

### Characteristics of the UL dataset from Koren et al (2024)

Sequencing was done in 2022 on R9.41, data comes from 8 separate Flowcells. Looking at the read length distribution the dataset available for download seems to be prefiltered since it does not contain shorter reads. However we found no information about the filtering tresholds and parameters.

### Comparison dataset overview

Overview of the first two flowcells we sequenced. This is unprocessed reads for our two flowcells and data "as downloaded" from Koren et al 2024.


| filename     | yield        | N50   | max readlen | yield>100kb | yield>200kb | reads > 1Mb |
| :----------- | :----------- | :---- | :----------- | :---------- | :---------- | :---------- |
| published.UL | 496106934638 | 81748 | 1.9 Mb       | 190 Gb      | 31 Gb       | 45          |
| run_04399.UL | 113914796731 | 67652 | 1.7 Mb       | 38 Gb       | 10 Gb       | 29          |
| run_04400.UL | 107609849976 | 77925 | 2.0 Mb       | 40 Gb       | 9 Gb        | 54          |


### Evaluation first UL sequencing runs.

Looking at the sequencing performance over time we see a high and steady output over the whole runtime. 
![Sequencing over time for first two UL runs](../doc/img/UL_sequencing_firstruns_yield.PNG)

Our sequencing yield is good, we can assume that it is higher then the published dataset where 8 Flowcells were combined to reach 500Gb after filtering leading to an average of 62.5Gb/Flowcell. We reached 230Gb with two flowcells, however before filtering.

### Read length and quality distribution

Before filtering:
![Read quality and length density](img/UL_comparison.readstats_density.png)

After filtering to 50x coverage (#todo):

Mean read quality is substantially better with the new R10 pores, the variance is broader with some low quality reads. This observations could stem from some pre filtering of the published dataset. 
Read length distribution of our UL reads is shifted slightly to slightly lower reads, reflected by 

There is no clear difference in read or quality length distribution between our two runs, the extraction methods seem to perform very similar. 

### Estimation of flowcells needed per sample for T2T-Assembly
