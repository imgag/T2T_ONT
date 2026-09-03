# Duplex basecalling

Duplex basecalling was evaluated early in the project as an alternative to HERRO correction. HERRO correction was ultimately selected for this study as it produces higher quality corrected reads with fewer flowcells required.

For reference, the standalone duplex pipeline is available at: https://github.com/imgag/NCCT_scripts

The duplex pipeline is also integrated into the T2T-ONT Snakemake workflow (`workflow/rules/basecalling.smk`).

## Duplex rates

Duplex rates of approximately 25% were observed on LSK-109 flowcells without specific optimizations for duplex sequencing.
