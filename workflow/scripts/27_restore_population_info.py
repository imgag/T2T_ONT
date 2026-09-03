#!/usr/bin/env python3
"""
Restore population information to merged PSAM file
"""
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Restore population info to merged PSAM')
    parser.add_argument('--merged-psam', required=True, help='Merged PSAM file (missing pop info)')
    parser.add_argument('--ref-psam', required=True, help='Reference PSAM file (with pop info)')
    parser.add_argument('--sample-psam', required=True, help='Sample PSAM file')
    parser.add_argument('--output', required=True, help='Output enhanced PSAM file')
    
    args = parser.parse_args()
    
    print("Reading merged PSAM file...")
    merged_df = pd.read_csv(args.merged_psam, sep='\t')
    
    print("Reading reference PSAM file...")
    ref_df = pd.read_csv(args.ref_psam, sep='\t')
    
    print("Reading sample PSAM file...")
    sample_df = pd.read_csv(args.sample_psam, sep='\t')
    
    # Create a lookup dictionary for population information
    pop_lookup = {}
    
    # Add reference population info
    if 'SUPERPOP' in ref_df.columns and 'POP' in ref_df.columns:
        for _, row in ref_df.iterrows():
            iid = row['IID']
            pop_lookup[iid] = {
                'POP': row['POP'],
                'SUPERPOP': row['SUPERPOP']
            }
    
    # Add sample population info (if available)
    if 'SUPERPOP' in sample_df.columns and 'POP' in sample_df.columns:
        for _, row in sample_df.iterrows():
            iid = row['IID']
            pop_lookup[iid] = {
                'POP': row['POP'],
                'SUPERPOP': row['SUPERPOP']
            }
    elif 'POP' in sample_df.columns:
        for _, row in sample_df.iterrows():
            iid = row['IID']
            pop_lookup[iid] = {
                'POP': row['POP'],
                'SUPERPOP': 'QUERY'  # Mark query samples
            }
    else:
        # No population info in sample file - mark all as query
        for _, row in sample_df.iterrows():
            iid = row['IID']
            pop_lookup[iid] = {
                'POP': 'QUERY',
                'SUPERPOP': 'QUERY'
            }
    
    # Add population columns to merged data
    merged_df['POP'] = merged_df['IID'].map(lambda x: pop_lookup.get(x, {}).get('POP', 'UNK'))
    merged_df['SUPERPOP'] = merged_df['IID'].map(lambda x: pop_lookup.get(x, {}).get('SUPERPOP', 'UNK'))
    
    # Write enhanced PSAM file
    print(f"Writing enhanced PSAM file with {len(merged_df)} samples...")
    merged_df.to_csv(args.output, sep='\t', index=False)
    
    # Report statistics
    print(f"Population distribution:")
    pop_counts = merged_df['SUPERPOP'].value_counts()
    for pop, count in pop_counts.items():
        print(f"  {pop}: {count} samples")
    
    print(f"Enhanced PSAM file created: {args.output}")

if __name__ == '__main__':
    main()