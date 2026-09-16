#!/usr/bin/env python3
"""
Collect and merge modkit_summary files into a single table.
Extracts sample name, reference type, and haplotype information.
"""

import sys
import os
import re
from pathlib import Path
import pandas as pd

def parse_modkit_summary(filepath):
    """Parse a single modkit_summary.tsv file."""
    data = {}
    modifications = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Parse header metadata
            if line.startswith('# bases'):
                data['base'] = line.split()[2]
            elif line.startswith('# total_reads_used'):
                data['total_reads_used'] = int(line.split()[2])
            elif line.startswith('# count_reads_'):
                data['count_reads'] = int(line.split()[2])
            elif line.startswith('# pass_threshold_'):
                data['pass_threshold'] = float(line.split()[2])
            
            # Parse modification data (skip header line with "base code pass_count...")
            elif not line.startswith('#') and not line.startswith('base'):
                parts = line.split()
                if len(parts) == 6:
                    mod_data = {
                        'base': parts[0],
                        'code': parts[1],
                        'pass_count': int(parts[2]),
                        'pass_frac': float(parts[3]),
                        'all_count': int(parts[4]),
                        'all_frac': float(parts[5])
                    }
                    modifications.append(mod_data)
    
    return data, modifications

def extract_info_from_filename(filepath):
    """Extract sample, ref, and haplotype from filename."""
    filename = Path(filepath).stem
    # Pattern: {sample}.{ref}.modkit_summary or {sample}.{ref}.hp{N}.modkit_summary
    
    parts = filename.split('.')
    sample = parts[0]
    
    # Check if it's a haplotype-specific file
    if 'hp1' in filename or 'hp2' in filename:
        # Pattern: sample.ref.hp1 or sample.ref.hp2
        haplotype = 'hp1' if 'hp1' in filename else 'hp2'
        ref = 'ref'
    else:
        # Pattern: sample.ref or sample.asm
        ref = parts[1]
        haplotype = 'combined'
    
    return sample, ref, haplotype

def main(input_files, output_file):
    """Merge all modkit_summary files into a single table."""
    
    all_data = []
    
    for filepath in input_files:
        if not os.path.exists(filepath):
            print(f"Warning: File not found: {filepath}", file=sys.stderr)
            continue
        
        sample, ref, haplotype = extract_info_from_filename(filepath)
        metadata, modifications = parse_modkit_summary(filepath)
        
        # Create a row for each modification type
        for mod in modifications:
            row = {
                'sample': sample,
                'reference': ref,
                'haplotype': haplotype,
                'total_reads_used': metadata.get('total_reads_used', 'NA'),
                'count_reads': metadata.get('count_reads', 'NA'),
                'pass_threshold': metadata.get('pass_threshold', 'NA'),
                'base': mod['base'],
                'mod_code': mod['code'],
                'pass_count': mod['pass_count'],
                'pass_frac': mod['pass_frac'],
                'all_count': mod['all_count'],
                'all_frac': mod['all_frac']
            }
            all_data.append(row)
    
    # Create DataFrame and save
    df = pd.DataFrame(all_data)
    
    # Sort by sample, reference, haplotype, then modification code
    df = df.sort_values(['sample', 'reference', 'haplotype', 'mod_code'])
    
    # Save to TSV
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Successfully merged {len(input_files)} files into {output_file}")
    print(f"Total rows: {len(df)}")
    print(f"Samples: {df['sample'].nunique()}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python collect_methylation_summary.py <output_file> <input_file1> <input_file2> ...", file=sys.stderr)
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    main(input_files, output_file)