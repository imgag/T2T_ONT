mamba activate activate /mnt/storage2/users/ahgrosc1/environments/envs/snakemake/envs/repeatmasker/

RepeatMasker \
    -species human \
    -dir analysis_other/repeatmasker/TUE_02_03UL \
    -pa 64 \
    -gff \
    -html \
    assembly/output/verkko/TUE_02_03UL/assembly.fasta

# Run the R script to analyze the RepeatMasker output
Rscript workflow/scripts/16_analyze_repeatmasker.R \
    --input analysis_other/repeatmasker/TUE_02_03UL/assembly.fasta.out \
    --output analysis_other/repeatmasker/rm_summary/TUE_02_03UL 

