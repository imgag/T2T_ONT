rule change_coverage:
    input:
        fq_dup = "data/combined/Duplex/1_3_23_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass_all.fastq",
        fq_ul = "data/combined/UL/03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup_all.fastq",
	params:
		target_cov_dup = config[target_coverage_hq]
		target_cov_ul = config[target_coverage_ul]
		target_base_dup = config[target_coverage_hq] * 3200000000
		target_base_ul = config[target_coverage_ul] * 3200000000
    output:
        out_dup = "data/combined/Duplex/1_3_23_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass_adj_cov.fastq.gz"
		out_ul = "data/combined/UL/03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup_adj_cov.fastq.gz"
    conda:
        "../env/filtlong.yml"
    log:
        "logs/filtlong.log"
    shell:
        """
        filtlong --min_length 10000 --target_bases {params.target_base_dup} --length_weight 10 1_3_23_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass_all.fastq | gzip > 1_3_23_R1041_Duplex_HG002_Dorado_v0.1.1_400bps_sup_stereo_duplex_pass_adj_cov.fastq.gz
        filtlong --min_mean_q 7 --min_length 50000 --target_bases {params.target_base_ul} --mean_q_weight 10 03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup.fastq | gzip > 03_08_22_R941_HG002_Guppy_6.1.2_5mc_cg_prom_sup_adj_cov.fastq.gz
        """
