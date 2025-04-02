# Run Merqury on polished and unpolished assemblies (not phased) against short reads

# Activate conda env
conda activate /mnt/storage2/users/ahgrosc1/environments/envs/snakemake/envs/merqury
export PATH=$PATH:"$CONDA_PREFIX"/share/merqury/eval

# Convert cram to fastqs
samtools fastq /mnt/storage2/projects/diagnostic/Genome_Diagnostik/Sample_DX203429_02/DX203429_02.cram \
    | pigz -p 12 > data/short_read/DX203429_02.fastq.gz

# Build k-mer db from Illumnina data
meryl k=21 count data/short_read/DX203429_02.fastq.gz output analysis_other/merqury_shortread/DX203429_02.meryl

# Run Merqury
INPUT_MERYL=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/DX203429_02.meryl
INPUT_PAT_FA=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/output/verkko_unphased/TUE_02/assembly.fasta
OUTPUT_PREFIX=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/TUE_02/merqury
LOG_FILE=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/TUE_02/merqury.log
pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
qv.sh \
    $INPUT_MERYL \
    $INPUT_PAT_FA \
    $OUTPUT_PREFIX \
    >> $LOG_FILE 2>&1
popd >> $LOG_FILE 2>&1


INPUT_MERYL=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/DX203429_02.meryl
INPUT_PAT_FA=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/assembly/output/verkko_unphased/TUE_02/assembly_polished.fasta
OUTPUT_PREFIX=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/TUE_02/merqury_polished
LOG_FILE=/mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/analysis_other/merqury_shortread/TUE_02/merqury_polished.log
pushd $(dirname $OUTPUT_PREFIX) >$LOG_FILE 2>&1
qv.sh \
    $INPUT_MERYL \
    $INPUT_PAT_FA \
    $OUTPUT_PREFIX \
    >> $LOG_FILE 2>&1
popd >> $LOG_FILE 2>&1
