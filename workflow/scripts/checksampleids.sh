#!/bin/bash
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/checksampleids.sh

# Extract sample IDs from VCF and PSAM files
echo "Extracting sample IDs..."

# Get sample IDs from VCF (header line with sample names) - PRESERVE ORDER
bcftools query -l analysis_other/ancestry/processed_vcf/1000G/1000G_merged_T2T.vcf.gz > vcf_samples_ordered.txt
echo "VCF samples: $(wc -l < vcf_samples_ordered.txt)"

# Get sample IDs from PSAM (column 2, skip header) - PRESERVE ORDER
tail -n +2 analysis_other/ancestry/plink/reference/1000G_phase3_T2T.psam | cut -f2 > psam_samples_ordered.txt
echo "PSAM samples: $(wc -l < psam_samples_ordered.txt)"

# Create sorted versions for set operations
sort vcf_samples_ordered.txt > vcf_samples_sorted.txt
sort psam_samples_ordered.txt > psam_samples_sorted.txt

echo "=== SET COMPARISON ==="

# Find samples only in VCF
comm -23 vcf_samples_sorted.txt psam_samples_sorted.txt > only_in_vcf.txt
echo "Only in VCF: $(wc -l < only_in_vcf.txt)"
if [ -s only_in_vcf.txt ]; then
    echo "First 10 samples only in VCF:"
    head -10 only_in_vcf.txt
fi

# Find samples only in PSAM
comm -13 vcf_samples_sorted.txt psam_samples_sorted.txt > only_in_psam.txt
echo "Only in PSAM: $(wc -l < only_in_psam.txt)"
if [ -s only_in_psam.txt ]; then
    echo "First 10 samples only in PSAM:"
    head -10 only_in_psam.txt
fi

# Find common samples
comm -12 vcf_samples_sorted.txt psam_samples_sorted.txt > common_samples.txt
echo "Common samples: $(wc -l < common_samples.txt)"

echo ""
echo "=== ORDER COMPARISON ==="

# Check if the order is exactly the same
if diff -q vcf_samples_ordered.txt psam_samples_ordered.txt > /dev/null; then
    echo "✅ Sample order is IDENTICAL between VCF and PSAM"
else
    echo "❌ Sample order is DIFFERENT between VCF and PSAM"
    
    # Show first few differences
    echo ""
    echo "First 10 differences in order:"
    echo "Position | VCF Sample | PSAM Sample"
    echo "---------|------------|------------"
    
    # Compare line by line and show first 10 differences
    paste -d'|' <(nl vcf_samples_ordered.txt) <(nl psam_samples_ordered.txt) | \
    awk -F'|' '$2 != $4 { print $1 " | " $2 " | " $4; count++; if (count >= 10) exit }'
    
    echo ""
    echo "To see all differences, run:"
    echo "diff -u vcf_samples_ordered.txt psam_samples_ordered.txt | head -50"
fi

echo ""
echo "=== DETAILED STATISTICS ==="

# Check if both files have the same number of samples
VCF_COUNT=$(wc -l < vcf_samples_ordered.txt)
PSAM_COUNT=$(wc -l < psam_samples_ordered.txt)

if [ "$VCF_COUNT" -eq "$PSAM_COUNT" ]; then
    echo "✅ Same number of samples: $VCF_COUNT"
else
    echo "❌ Different number of samples: VCF=$VCF_COUNT, PSAM=$PSAM_COUNT"
fi

# Check for duplicates
VCF_UNIQUE=$(sort vcf_samples_ordered.txt | uniq | wc -l)
PSAM_UNIQUE=$(sort psam_samples_ordered.txt | uniq | wc -l)

echo "Unique samples - VCF: $VCF_UNIQUE, PSAM: $PSAM_UNIQUE"

if [ "$VCF_COUNT" -ne "$VCF_UNIQUE" ]; then
    echo "❌ VCF has duplicate samples!"
    echo "Duplicates in VCF:"
    sort vcf_samples_ordered.txt | uniq -d | head -5
fi

if [ "$PSAM_COUNT" -ne "$PSAM_UNIQUE" ]; then
    echo "❌ PSAM has duplicate samples!"
    echo "Duplicates in PSAM:"
    sort psam_samples_ordered.txt | uniq -d | head -5
fi

echo ""
echo "=== SAMPLE PREVIEW ==="
echo "First 5 VCF samples:"
head -5 vcf_samples_ordered.txt

echo ""
echo "First 5 PSAM samples:"
head -5 psam_samples_ordered.txt

echo ""
echo "Last 5 VCF samples:"
tail -5 vcf_samples_ordered.txt

echo ""
echo "Last 5 PSAM samples:"
tail -5 psam_samples_ordered.txt

# Cleanup
rm vcf_samples_ordered.txt psam_samples_ordered.txt vcf_samples_sorted.txt psam_samples_sorted.txt only_in_vcf.txt only_in_psam.txt common_samples.txt

echo ""
echo "Sample order check completed!"