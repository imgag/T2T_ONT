# Tasks

- [ ] Check how a combination of Duplex and HERRO reads performs?
    - [ ] Is the coverage evenly distributed?
    - [ ] Is it possible to exclude reads that were used in Duplex generation from HERRO calling
    - [ ] What is the max benefitial addition of HERRO reads?
- [ ] Fix the basecalling to include methylation (!) in the real samples. Only needed for Ultra-Long. 
- [ ] Compare assembly results with publication
- [ ] Generate PoreC statistics with the “evaluate_contacts” executable provided in the GFAse package
- [ ] Generate Read and Mapping stats for PoreC
- [ ] Documentation
    - [ ] Generate UL Read QC presentation
    - [ ] Write methods section for assembly approach
    - [ ] Write overview for used published datasets
- [ ] Create Beautiful figures for publication or print
    - Stained Glass [Link](https://resgen.io/paper-data/T2T/views/MtjcVgrlQmymnHIvdck5-g)
    - Coloured GFAse assembly plot output
    - Karyotype plots with annotation of telomeres, centromeres, phasing blocks, synteny
- [ ] Investigate whats going on with HERRO read lenghts and how it is affectign subsampling. This plot looks strange: `doc/img/HQ_herro_published_subsampling.read_stats.png`
- [ ] Check if we used the same data for HERRO as in the paper. Data depended errors in the model? Can we achieve results from HERRO publication?
- [ ] Check if filtlong subsampling with length preference is a problem for subsampling HERRO and DUPLEX reads in the benchmark

- [X] Make PORCE optional, dont run scaffolding without it
- [X] Telomere detection (seqtk telo ?, take from verrko results, )
- [X] Is our sample TUE_01 the same HG002 from published dataset? No its a personal genome from Germany

- [X] Use [Fastcat](https://github.com/epi2me-labs/fastcat) bamqc for bamstats
- [X] Investigate if what is inside the [T2T bundle[(https://labs.epi2me.io/lc2024_t2 )] from ONT. Maybe includes pipeline?
- [X] Include [asmgene](https://github.com/lh3/minimap2/issues/355) QC for gene completenes stats
- [X] Generate sequencing stats on raw data files
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
    - [X] Compare alignment QC with Duplex
- [X] Setup VERRKO pipeline
- [X] Run VERRKO assembly
    - [X] Duplex reads
    - [X] HERRO corrected reads
    - [X] Normal/no duplex reads
