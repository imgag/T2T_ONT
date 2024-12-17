---
marp: true
paginate: true
theme: gaia
title: 2024-12-17 IMAGAG longread meeting
--- 

# IMGAG longread meeting

17.12.2024
1st status update T2T
Caspar and Prithivi


---
### Important Resources
- Koren et al. 2024: Gapless assembly of complete human and plant chromosomes using only nanopore sequencing (**Verkko**)
- Šikić et al. 2024: Telomere-to-telomere phased genome assembly using error-corrected Simplex nanopore reads (**Herro**)
- Li and Durbin 2024: Genome assembly in the telomere-to-telomere era (**Review**)
- ONT T2T Assembly Blogpost : https://community.nanoporetech.com/t/ont-t2t-assembly-blogpost/10166

--- 

### Nanopore Read Types
- **Simplex**: Default ONT reads, LSK Library (R10.3)
- **Duplex**: Simplex reads, but with duplex chemistry (R10.3)
- **Herro**: Simplex or UL reads, corrected with ML algorithm (R10.3)
- **UL**: Special Library Prep, focus on max read length (R10.4)
- **PoreC**: 
- **Polishing Kit**: Simplex reads, but with Polishing Kit chemistry 

--- 

### Public Datasets

Giab Sample **HG002** From Koren et al (2024): 
- 15 FC Simplex R10
- Duplex reads called from simplex reads 
- 8 FC R9.4 UltraLong
- PoreC

--- 

### Published Pipelines
- ONT Recommended (Verkko with Herro)
![bg right:60% h:500](../img/asm-workflow-ont-blog.png)
- Verkko with Duplex (Koren et al. 2024)

--- 

- Datasets
    - Published
    - First Seqs
--- 

### Benchmarking with published data

Run on reduced dataset: Filtered on chr19

![Overview chr19 in CHM13](../img/chr19_T2T_ref.png)

- *"Easy"* chromosome
- Fast computation

--- 

### Comparison table
```
 LR↓/UL→ |         50x         |      70x     |     90x      |
         +---------------------+--------------+--------------+
     15x |     duplex          |       --     |      --      |
     20x |     duplex          |       --     |      --      |
     25x |     duplex          |       --     |      --      |
     35x | 15x DU + 20x herro  |       --     |      --      |
     35x | duplex/herro        | duplex/herro | duplex/herro |
     60x |     herro           |     herro    |    herro     |
    120x |     herro           |     herro    |    herro     |
```
---


1. Filtering on fastq. Do we reach target coverage? 

![bg contain](../img/subsampling_coverages.cov_boxplot.png)


---
 
## Impact UL Coverage


---

## Impact HQ Coverage


---

35x Herro 50x UL chr19 <-----------------> 35x Duplex 50x UL chr19

![bg h:550](../../assembly/qc/phased_verkko/published_chr19_herro_35x_UL50x/bandage_graph.png)
![bg h:550](../../assembly/qc/phased_verkko/published_chr19_duplex_UL50x/bandage_graph.png)

---

### Completed T2T ?
no postprocessing, T2T detection with script, looking for TTAGGG

Duplex, UL70x, Verkko
Haplotype1: **16** chromosomes T2T
Haplotype2: **14** chromosomes T2T

Duplex, UL70x, Gfase
Haplotype1: **8** chromosomes T2T
Haplotype2: **0** chromosomes T2T

--- 

### HQ Read comparison
    - Scaffolder comparison


--- 

Our data: UL

![bg h:680](../img/UL_comparison_50x.read_stats.png)

---
### Our data: Duplex



---
## 5. Discussions
    - T2T requirement discussions
        - Genome completeness
        - Planned number of Flowcells 
    - Lab Establishment Priorities
        - UL
        - Duplex
        - PoreC
        - Polishing Kit
    - Bioinformatics Implementations