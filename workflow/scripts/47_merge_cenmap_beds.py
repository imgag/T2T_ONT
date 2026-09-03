#!/usr/bin/env python3
"""
Merge cenMAP BED files and add sample, chromosome, and haplotype columns.
"""

import argparse
import pandas as pd
import re
import sys


def parse_name_column(name_str):
    """
    Parse the name column (e.g., 'T2T00_1_chr10_haplotype1-0000008' or 'T2T00_1_rc-chr13_haplotype1-0000016')
    to extract sample, chromosome, and haplotype.
    
    Returns: (sample, chromosome, haplotype)
    """
    # Pattern to match: sample_chr(X)_haplotype or sample_rc-chr(X)_haplotype
    # Sample can be T2T00_1 or T2T01 (with optional underscore and number)
    pattern = r'^(T2T\d+(?:_\d+)?)_(rc-)?chr([0-9XY]+)_(haplotype\d+)-\d+$'
    
    match = re.match(pattern, name_str)
    if match:
        sample = match.group(1)
        # rc prefix is ignored as strand already represents orientation
        chrom = f"chr{match.group(3)}"
        haplotype = match.group(4)
        return sample, chrom, haplotype
    else:
        # Fallback: try to extract what we can
        print(f"Warning: Could not parse name column: {name_str}", file=sys.stderr)
        return None, None, None


def process_bed_file(filepath):
    """
    Read and process a single BED file.
    """
    # BED files have no header
    # Columns: contig, start, end, name, score, strand, thickStart, thickEnd, itemRgb
    df = pd.read_csv(filepath, sep='\t', header=None,
                     names=['contig', 'start', 'end', 'name', 'score', 'strand', 
                            'thick_start', 'thick_end', 'item_rgb'])
    
    # Extract sample, chromosome, haplotype from name column
    parsed = df['name'].apply(parse_name_column)
    df['sample'] = parsed.apply(lambda x: x[0])
    df['chromosome'] = parsed.apply(lambda x: x[1])
    df['haplotype'] = parsed.apply(lambda x: x[2])
    
    # Calculate centromere length
    df['length'] = df['end'] - df['start']
    
    # Select and reorder columns (remove non-informative ones)
    # Keep: sample, chromosome, haplotype, contig, start, end, length, strand
    df_out = df[['sample', 'chromosome', 'haplotype', 'contig', 'start', 'end', 'length', 'strand']]
    
    return df_out


def main():
    parser = argparse.ArgumentParser(description='Merge cenMAP BED files')
    parser.add_argument('--input', '-i', nargs='+', required=True,
                        help='Input BED files')
    parser.add_argument('--output', '-o', required=True,
                        help='Output merged TSV file')
    
    args = parser.parse_args()
    
    all_data = []
    
    for filepath in args.input:
        try:
            df = process_bed_file(filepath)
            all_data.append(df)
        except Exception as e:
            print(f"Error processing {filepath}: {e}", file=sys.stderr)
            continue
    
    if all_data:
        combined = pd.concat(all_data, ignore_index=True)
        # Sort by sample, chromosome, haplotype
        combined = combined.sort_values(['sample', 'chromosome', 'haplotype'])
        combined.to_csv(args.output, sep='\t', index=False)
        print(f"Merged {len(args.input)} files into {args.output}", file=sys.stderr)
    else:
        print("No data to merge", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
