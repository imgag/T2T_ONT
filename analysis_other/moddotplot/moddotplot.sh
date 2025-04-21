# Reference dataset must be copied into the environment


HP1="assembly/output/verkko/TUE_02_03UL/assembly.haplotype1.fasta"
HP2="assembly/output/verkko/TUE_02_03UL/assembly.haplotype2.fasta"

moddotplot interactive -f analysis_other/moddotplot/TUE_02_03UL.haplotype1.chr13.fasta
ssh -N -f -L 8050:127.0.0.1:8050 srv023
