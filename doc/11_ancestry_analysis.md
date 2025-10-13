## Typical SNP Numbers in Ancestry Analysis

### **1. Global Ancestry Analysis**

| **Method**    | **Typical SNPs** | **Range** | **Examples from Literature**                                |
| ------------- | ---------------- | --------- | ----------------------------------------------------------- |
| **ADMIXTURE** | 100K-1M          | 50K-2M    | • 1000 Genomes: ~650K SNPs UK Biobank: ~700K SNPs      |
| **STRUCTURE** | 50K-500K         | 10K-1M    | • HapMap studies: ~300K SNPs • Population studies: 100K-500K |
| **iAdmix**    | 50K-300K         | 20K-500K  | • Original paper: ~200K SNPs • Your analysis: 27K SNPs   |
| **PCA-based** | 100K-1M          | 50K-2M    | • EIGENSTRAT: ~500K SNPs                                    |

### **2. Local Ancestry Analysis**

| **Method** | **Typical SNPs** | **Range** | **Notes**                         |
| ---------- | ---------------- | --------- | --------------------------------- |
| **RFMix**  | 500K-5M          | 100K-10M  | Needs dense coverage              |
| **LAMP**   | 300K-2M          | 100K-5M   | Sliding window approach           |
| **Gnomix** | 1M-10M           | 500K-20M  | Deep learning, more SNPs = better |



Literature Examples
## Landmark Studies:

1) 1000 Genomes Project (2015)
    - SNPs: ~650,000 high-quality SNPs
    - Used for: Global ancestry inference
    - Coverage: Genome-wide

2) Bryc et al. (2015) - African Americans
    - SNPs: ~750,000 SNPs
    - Tool: ADMIXTURE
    - Result: European ~73%, African ~24%

3) Moreno-Estrada et al. (2013) - Latin Americans
    - SNPs: ~300,000 SNPs
    - Tools: ADMIXTURE, RFMix
    - Multi-population analysis

4) Hellenthal et al. (2014) - Global mixing
    - SNPs: ~650,000 SNPs
    - Method: ChromoPainter/GLOBETROTTER
    - 95 populations worldwide

Recent High-Density Studies:

5) UK Biobank (2018-present)
   - SNPs: ~700,000 genotyped + imputed to millions

6) All of Us Research Program (2020-present)
   - SNPs: ~1M+ genotyped SNPs
   - Diverse populations


## Our SNP numbers

>`wc -l analysis_other/ancestry/global/iadmix/freq/reference_frequencies.txt` 
> 7397450


