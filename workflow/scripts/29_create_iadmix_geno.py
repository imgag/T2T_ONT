#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create iAdmix genotype file for a single sample from PLINK traw format')
    parser.add_argument('--plink-raw', required=True, help='PLINK .traw file')
    parser.add_argument('--sample', required=True, help='Sample ID to extract')
    parser.add_argument('--output', required=True, help='Output genotype file for iAdmix')
    
    args = parser.parse_args()
    
    print(f"Reading PLINK traw file for sample: {args.sample}")
    # Read PLINK traw file (transposed raw format)
    # Format: CHR SNP (C)M POS COUNTED ALT sample1 sample2 ...
    traw_df = pd.read_csv(args.plink_raw, sep='\t')
    
    print(f"Found {len(traw_df)} variants in traw file")
    
    # Get sample columns (everything after the first 6 columns)
    variant_info_cols = ['CHR', 'SNP', 'CM', 'POS', 'COUNTED', 'ALT']
    sample_cols = [col for col in traw_df.columns if col not in variant_info_cols]
    
    print(f"Available samples: {sample_cols}")
    
    # Check if requested sample exists
    if args.sample not in sample_cols:
        print(f"Error: Sample {args.sample} not found in traw file")
        print(f"Available samples: {sample_cols}")
        return 1
    
    print(f"Processing sample: {args.sample}")
    
    # Prepare output data
    output_data = []
    
    print("Converting genotype format...")
    
    # Process each variant
    for idx, row in traw_df.iterrows():
        if idx % 10000 == 0:
            print(f"Processed {idx} variants...")
        
        rsid = row['SNP']
        
        # Get genotype for this sample
        genotype_code = row[args.sample]
        
        # Convert PLINK coding (0, 1, 2, NA) to iAdmix format (AA, AB, BB)
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
        
        # Create row: rsid followed by single genotype
        output_data.append([rsid, iadmix_geno])
    
    # Write output file
    print(f"Writing genotype file to {args.output}")
    
    with open(args.output, 'w') as f:
        # Write header (rsid and sample name)
        f.write(f'rsid\t{args.sample}\n')
        
        # Write genotype data
        for row in output_data:
            f.write('\t'.join(row) + '\n')
    
    print("iAdmix genotype file created successfully!")
    print(f"Total variants: {len(output_data)}")
    print(f"Sample: {args.sample}")
    
    # Report genotype distribution
    genotypes = [row[1] for row in output_data]
    geno_counts = pd.Series(genotypes).value_counts()
    print("Genotype distribution:")
    for geno, count in geno_counts.items():
        print(f"  {geno}: {count}")

if __name__ == '__main__':
    main()