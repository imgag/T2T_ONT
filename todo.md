# Tasks

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
