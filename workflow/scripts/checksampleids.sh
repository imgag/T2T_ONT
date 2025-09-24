#!/bin/bash

# Extract sample IDs from VCF and PSAM files
echo "Extracting sample IDs..."

# Get sample IDs from VCF (header line with sample names)
bcftools query -l analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz | sort > vcf_samples.txt
echo "VCF samples: $(wc -l < vcf_samples.txt)"

# Get sample IDs from PSAM (column 2, skip header)
tail -n +2 analysis_other/ancestry/plink/reference/1000G_phase3_T2T.psam | cut -f2 | sort > psam_samples.txt
echo "PSAM samples: $(wc -l < psam_samples.txt)"

# Find samples only in VCF
comm -23 vcf_samples.txt psam_samples.txt > only_in_vcf.txt
echo "Only in VCF: $(wc -l < only_in_vcf.txt)"
if [ -s only_in_vcf.txt ]; then
    echo "First 10 samples only in VCF:"
    head -10 only_in_vcf.txt
fi

# Find samples only in PSAM
comm -13 vcf_samples.txt psam_samples.txt > only_in_psam.txt
echo "Only in PSAM: $(wc -l < only_in_psam.txt)"
if [ -s only_in_psam.txt ]; then
    echo "First 10 samples only in PSAM:"
    head -10 only_in_psam.txt
fi

# Find common samples
comm -12 vcf_samples.txt psam_samples.txt > common_samples.txt
echo "Common samples: $(wc -l < common_samples.txt)"

# Cleanup
rm vcf_samples.txt psam_samples.txt only_in_vcf.txt only_in_psam.txt common_samples.txt