mamba activate activate /mnt/storage2/users/ahgrosc1/environments/envs/snakemake/envs/repeatmasker/

RepeatMasker \
    -species human \
    -dir analysis_other/repeatmasker/TUE_02_03UL \
    -pa 64 \
    -gff \
    -html \
    assembly/output/verkko/TUE_02_03UL/assembly.fasta

awk -F'\t' 'BEGIN{OFS="\t"} $1!~/#/ {print $1,$4-1,$5,$3,$7}' analysis_other/repeatmasker/TUE_02_03UL/assembly.fasta.out.gff > analysis_other/repeatmasker/TUE_02_03UL/assembly.fasta.out.bed