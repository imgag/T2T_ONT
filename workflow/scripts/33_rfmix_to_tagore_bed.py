#!/usr/bin/env python3
"""
Convert RFMix MSP output to BED format for Tagore visualization.
Creates a single BED file with both haplotypes.
"""

import pandas as pd
import argparse
import os
import sys

# Population code mapping
POP_CODES = {
    '0': 'AFR',
    '1': 'AMR', 
    '2': 'EAS',
    '3': 'EUR',
    '4': 'SAS'
}

# Color scheme for populations (hex colors for Tagore)
POP_COLORS = {
    'AFR': '#FF6B6B',  # Red
    'AMR': '#4ECDC4',  # Teal
    'EAS': '#FFE66D',  # Yellow
    'EUR': '#95E1D3',  # Mint
    'SAS': '#C7CEEA'   # Lavender
}

def read_msp_file(msp_file):
    """
    Read RFMix MSP file, handling the header properly.
    Returns a DataFrame with proper column names.
    """
    with open(msp_file, 'r') as f:
        lines = f.readlines()
    
    # Find the column header line
    header_line = None
    data_start_idx = 0
    
    for i, line in enumerate(lines):
        if line.startswith('#chm'):
            header_line = line.strip()
            data_start_idx = i + 1
            break
    
    if not header_line:
        raise ValueError("Could not find column header line starting with #chm")
    
    # Parse column names and remove # from first column
    col_names = header_line.split('\t')
    col_names[0] = col_names[0].lstrip('#')
    
    # Read the data
    df = pd.read_csv(msp_file, sep='\t', skiprows=data_start_idx, 
                     names=col_names, header=None)
    
    return df

def convert_msp_to_tagore_bed(msp_file, sample_name, output_file):
    """
    Convert RFMix MSP output to Tagore BED format for a specific sample.
    Creates a single BED file with both haplotypes.
    
    Tagore BED format:
    #chr  start  stop  feature  size  color  chrCopy
    - feature: 0=rectangle, 1=circle, 2=triangle, 3=line
    - size: 0-1 (width of feature)
    - chrCopy: 1 or 2 (which haplotype)
    """
    
    # Read MSP file with proper header handling
    df = read_msp_file(msp_file)
    
    print(f"Loaded MSP file with {len(df)} segments", file=sys.stderr)
    print(f"Available columns (first 10): {df.columns.tolist()[:10]}", file=sys.stderr)
    
    # Get column names for this sample's haplotypes
    hap0_col = f"{sample_name}.0"
    hap1_col = f"{sample_name}.1"
    
    # Check if columns exist
    if hap0_col not in df.columns or hap1_col not in df.columns:
        # Try to find the sample in the columns
        sample_cols = [col for col in df.columns if sample_name in col]
        print(f"Sample columns found: {sample_cols}", file=sys.stderr)
        raise ValueError(f"Sample {sample_name} not found in MSP file. Columns with sample name: {sample_cols}")
    
    print(f"Processing sample {sample_name} with columns: {hap0_col}, {hap1_col}", file=sys.stderr)
    
    # Create BED entries for both haplotypes
    bed_rows = []
    
    for _, row in df.iterrows():
        chrom = row['chm']
        start = int(row['spos'])
        end = int(row['epos'])
        
        # Process both haplotypes
        for hap_idx, hap_col in enumerate([hap0_col, hap1_col], start=1):
            ancestry_code = str(int(row[hap_col]))
            ancestry = POP_CODES.get(ancestry_code, 'UNK')
            color = POP_COLORS.get(ancestry, '#CCCCCC')
            
            # Tagore BED format: chr, start, stop, feature, size, color, chrCopy
            bed_rows.append({
                '#chr': chrom,       # Note: using #chr as column name
                'start': start,
                'stop': end,
                'feature': 0,  # 0 = rectangle
                'size': 1,     # Full width
                'color': color,
                'chrCopy': hap_idx  # 1 or 2 for haplotype
            })
    
    # Create DataFrame and save
    bed_df = pd.DataFrame(bed_rows)
    
    # Sort by chromosome, position, and haplotype
    bed_df['chr_num'] = (bed_df['#chr']
                         .str.replace('chr', '', regex=False)
                         .replace('X', '23')
                         .replace('Y', '24')
                         .astype(int))
    bed_df = bed_df.sort_values(['chr_num', 'start', 'chrCopy']).drop('chr_num', axis=1)
    
    # Write with header - pandas will write #chr as the first column header
    bed_df.to_csv(output_file, sep='\t', index=False, header=True)
    
    print(f"Created {output_file} with {len(bed_df)} segments ({len(bed_df)//2} per haplotype)", 
          file=sys.stderr)
    
    # Print ancestry distribution
    ancestry_counts = {}
    for _, row in bed_df.iterrows():
        color = row['color']
        ancestry = [k for k, v in POP_COLORS.items() if v == color][0]
        ancestry_counts[ancestry] = ancestry_counts.get(ancestry, 0) + 1
    
    print(f"Ancestry distribution: {ancestry_counts}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Convert RFMix MSP output to Tagore BED format'
    )
    parser.add_argument('--msp', required=True, help='Input RFMix MSP file')
    parser.add_argument('--sample', required=True, help='Sample ID to extract')
    parser.add_argument('--output', required=True, 
                       help='Output BED file for Tagore')
    
    args = parser.parse_args()
    
    # Convert MSP to Tagore BED
    convert_msp_to_tagore_bed(args.msp, args.sample, args.output)

if __name__ == '__main__':
    main()