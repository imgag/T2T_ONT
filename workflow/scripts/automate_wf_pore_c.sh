cd bin
curl -s https://get.sdkman.io | bash
chmod +x nextflow
export PATH=$PATH:/mnt/storage2/users/ahthapp1/real_deal/github_projects/T2T_ONT/bin
nextflow run bin/wf-pore-c/main.nf --fastq /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/data/raw/24070/24070LRa_PoreC_04921/20250324_1531_1C_PBC81151_52756a5e/fastq_pass \ 
	--ref data/ref/T2T-CHM13.v2.fasta --out_dir /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/output-wf-pore-c --hi_c True --bed True --threads 19 -profile singularity --coverage True --chromunity True --pairs True --mcool True
