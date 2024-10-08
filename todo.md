# Tasks

- [ ] Prepare published test datasets:
    - [X] Download datasets (see Readme)
    - [X] Merge bamfiles
    - [ ] Convert to suitable input for Verrko (.fastq.gz)
- [ ] Create subsampled testdataset:
    - [ ] Map UL and Duplex against T2T assembly
    - [ ] Extract all reads matching chr19 (small, short telomeres and centromeres).
    ![chr1 from T2T ref assembly](doc/img/chr19_T2T_ref.png)
- [ ] Setup HERRO pipeline
    - [ ] Process SIMPLEX reads with HERRO
    - [ ] Compare alignment QC with Duplex
- [ ] Setup VERRKO pipeline
- [ ] Run VERRKO assembly
    - [ ] Duplex reads
    - [ ] HERRO corrected reads
    - [ ] Normal/no duplex reads
- [ ]  Compare assembly results with publication
