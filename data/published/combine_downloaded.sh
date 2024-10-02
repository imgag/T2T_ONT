find  downloads/Duplex -name '*_pass.bam' -exec samtools merge -o data/published/11_15_22_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass.bam {} +
find  downloads/UL -name '*.bam.gz' -exec gunzip {} \;
find  downloads/UL -name '*.bam' -exec samtools merge -o data/published/03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup.bam {} +
