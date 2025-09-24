#!/usr/bin/env python3
"""
Filter PSAM file to only include samples present in VCF file
"""
import argparse
import subprocess
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Filter PSAM to match VCF samples')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--metadata', required=True, help='Input metadata file (TSV)')
    parser.add_argument('--psam', required=True, help='Output PSAM file')
    
    args = parser.parse_args()
    
    # Get VCF sample IDs using bcftools
    print("Extracting sample IDs from VCF...", file=sys.stderr)
    try:
        result = subprocess.run(['bcftools', 'query', '-l', args.vcf], 
                               capture_output=True, text=True, check=True)
        vcf_samples = set(result.stdout.strip().split('\n'))
        print(f"Found {len(vcf_samples)} samples in VCF", file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Read metadata file
    print("Reading metadata file...", file=sys.stderr)
    try:
        metadata = pd.read_csv(args.metadata, sep='\t')
        print(f"Found {len(metadata)} samples in metadata", file=sys.stderr)
    except Exception as e:
        print(f"Error reading metadata: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Filter metadata to samples present in VCF
    # Assuming Individual ID is in column 1 (0-indexed)
    iid_column = metadata.columns[1]
    filtered_metadata = metadata[metadata[iid_column].isin(vcf_samples)]
    print(f"Found {len(filtered_metadata)} matching samples", file=sys.stderr)
    
    # Write PSAM file
    print("Writing PSAM file...", file=sys.stderr)
    with open(args.psam, 'w') as f:
        # Write PSAM header
        f.write('#FID\tIID\tPAT\tMAT\tSEX\tPHENO\n')
        
        # Write filtered data (first 6 columns)
        for _, row in filtered_metadata.iterrows():
            f.write(f'{row.iloc[0]}\t{row.iloc[1]}\t{row.iloc[2]}\t{row.iloc[3]}\t{row.iloc[4]}\t{row.iloc[5]}\n')
    
    # Report stats
    print(f"Original metadata samples: {len(metadata)}", file=sys.stderr)
    print(f"VCF samples: {len(vcf_samples)}", file=sys.stderr) 
    print(f"Filtered PSAM samples: {len(filtered_metadata)}", file=sys.stderr)
    
    # Report mismatches for debugging
    metadata_samples = set(metadata[iid_column])
    only_in_metadata = metadata_samples - vcf_samples
    only_in_vcf = vcf_samples - metadata_samples
    
    if only_in_metadata:
        print(f"Samples in metadata but not in VCF: {len(only_in_metadata)}", file=sys.stderr)
        if len(only_in_metadata) <= 10:
            print(f"  {list(only_in_metadata)}", file=sys.stderr)
    
    if only_in_vcf:
        print(f"Samples in VCF but not in metadata: {len(only_in_vcf)}", file=sys.stderr)
        if len(only_in_vcf) <= 10:
            print(f"  {list(only_in_vcf)}", file=sys.stderr)

if __name__ == '__main__':
    main()