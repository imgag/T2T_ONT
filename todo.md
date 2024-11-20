# Tasks

- [ ] Use [Fastcat](https://github.com/epi2me-labs/fastcat) bamqc for bamstats
- [ ] Investigate if what is inside the [T2T bundle[(https://labs.epi2me.io/lc2024_t2 )] from ONT. Maybe includes pipeline?
- [X] Include [asmgene](https://github.com/lh3/minimap2/issues/355) QC for gene completenes stats
- [ ] Generate sequencing stats on raw data files
- [X] Try different Ultralong coverages

- [X] Prepare published test datasets:
    - [X] Download datasets (see Readme)
    - [X] Merge bamfiles
    - [X] Convert to suitable input for Verrko (.fastq.gz)
- [X] Create subsampled testdataset:
    - [X] Map UL and Duplex against T2T assembly
    - [X] Extract all reads matching chr19 (small, short telomeres and centromeres).
    ![chr1 from T2T ref assembly](doc/img/chr19_T2T_ref.png)
- [X] Setup HERRO pipeline
    - [X] Process SIMPLEX reads with HERRO
    - [ ] Compare alignment QC with Duplex
- [X] Setup VERRKO pipeline
- [X] Run VERRKO assembly
    - [X] Duplex reads
    - [X] HERRO corrected reads
    - [X] Normal/no duplex reads
- [ ]  Compare assembly results with publication


Six sentence abstract
Meeting in Bremen
