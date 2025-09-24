#!/usr/bin/env python3
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/29_create_iadmix_geno.py

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create iAdmix genotype file from PLINK traw format')
    parser.add_argument('--plink-raw', required=True, help='PLINK .traw file')
    parser.add_argument('--output', required=True, help='Output genotype file for iAdmix')
    
    args = parser.parse_args()
    
    print("Reading PLINK traw file...")
    # Read PLINK traw file (transposed raw format)
    # Format: CHR SNP (C)M POS COUNTED ALT sample1 sample2 ...
    traw_df = pd.read_csv(args.plink_raw, sep='\t')
    
    print(f"Found {len(traw_df)} variants and {len(traw_df.columns) - 6} samples")
    
    # Get sample columns (everything after the first 6 columns)
    variant_info_cols = ['CHR', 'SNP', 'CM', 'POS', 'COUNTED', 'ALT']
    sample_cols = [col for col in traw_df.columns if col not in variant_info_cols]
    
    print(f"Sample columns: {sample_cols[:5]}..." if len(sample_cols) > 5 else f"Sample columns: {sample_cols}")
    
    # Prepare output data
    output_data = []
    
    print("Converting genotype format...")
    
    # Process each variant
    for idx, row in traw_df.iterrows():
        if idx % 1000 == 0:
            print(f"Processed {idx} variants...")
        
        rsid = row['SNP']
        counted_allele = row['COUNTED']
        alt_allele = row['ALT']
        
        # Create row for this variant
        variant_row = [rsid]
        
        # Process each sample
        for sample in sample_cols:
            genotype_code = row[sample]
            
            # Convert PLINK coding (0, 1, 2, NA) to iAdmix format (AA, AB, BB)
            # Simplified approach: use A/B notation since iAdmix mainly cares about the pattern
            if pd.isna(genotype_code) or genotype_code == -9:
                iadmix_geno = "NN"  # Missing genotype
            elif genotype_code == 0:
                iadmix_geno = "AA"  # Homozygous reference
            elif genotype_code == 1:
                iadmix_geno = "AB"  # Heterozygous
            elif genotype_code == 2:
                iadmix_geno = "BB"  # Homozygous alternate
            else:
                iadmix_geno = "NN"  # Unknown coding
            
            variant_row.append(iadmix_geno)
        
        output_data.append(variant_row)
    
    # Write output file
    print(f"Writing genotype file to {args.output}")
    
    with open(args.output, 'w') as f:
        # Write header
        header = ['rsid'] + sample_cols
        f.write('\t'.join(header) + '\n')
        
        # Write genotype data
        for row in output_data:
            f.write('\t'.join(row) + '\n')
    
    print("iAdmix genotype file created successfully!")
    print(f"Total variants: {len(output_data)}")
    print(f"Total samples: {len(sample_cols)}")

if __name__ == '__main__':
    main()