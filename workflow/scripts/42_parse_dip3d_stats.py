#!/usr/bin/env python3
"""
Parse dip3d statistics files and combine them into a single table.

This script extracts key metrics from dip3d haplotagging statistics files
(both imputed and snp-tagged variants) and combines them into a TSV table.
"""

import re
import sys
import argparse
import pandas as pd
from pathlib import Path


def parse_dip3d_stats(file_path, sample_name, tag_type):
    """
    Parse a single dip3d stats file and extract key metrics.
    
    Args:
        file_path: Path to the dip3d stats file
        sample_name: Name of the sample
        tag_type: Type of tagging (imputed or snp-tagged)
    
    Returns:
        Dictionary containing extracted statistics
    """
    stats = {
        'sample': sample_name,
        'tag_type': tag_type
    }
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract total reads
    match = re.search(r'Total reads: (\d+)', content)
    if match:
        stats['total_reads'] = int(match.group(1))
    
    # Extract reads with phased frags
    match = re.search(r'Reads with phased frags: (\d+) \(([0-9.]+)%\)', content)
    if match:
        stats['reads_with_phased_frags'] = int(match.group(1))
        stats['reads_with_phased_frags_pct'] = float(match.group(2))
    
    # Extract total frags
    match = re.search(r'Total frags: (\d+)', content)
    if match:
        stats['total_frags'] = int(match.group(1))
    
    # Extract phased frags
    match = re.search(r'Phased frags: (\d+) \(([0-9.]+)%\)', content)
    if match:
        stats['phased_frags'] = int(match.group(1))
        stats['phased_frags_pct'] = float(match.group(2))
    
    # Extract H1 frags
    match = re.search(r'H1 frags: (\d+) \(([0-9.]+)%\)', content)
    if match:
        stats['h1_frags'] = int(match.group(1))
        stats['h1_frags_pct'] = float(match.group(2))
    
    # Extract H2 frags
    match = re.search(r'H2 frags: (\d+) \(([0-9.]+)%\)', content)
    if match:
        stats['h2_frags'] = int(match.group(1))
        stats['h2_frags_pct'] = float(match.group(2))
    
    # Extract contact stats for different distance thresholds
    contact_sections = re.findall(
        r'Contact stats for contact dist <= (\d+)\s+Total contacts: (\d+)\s+Phased contacts: (\d+) \(([0-9.]+)%\)\s+H1 contacts: (\d+) \(([0-9.]+)%\)\s+H2 contacts: (\d+) \(([0-9.]+)%\)\s+H-trans contacts: (\d+) \(([0-9.]+)%\)',
        content
    )
    
    for section in contact_sections:
        dist, total, phased, phased_pct, h1, h1_pct, h2, h2_pct, htrans, htrans_pct = section
        prefix = f'contacts_dist_{dist}'
        stats[f'{prefix}_total'] = int(total)
        stats[f'{prefix}_phased'] = int(phased)
        stats[f'{prefix}_phased_pct'] = float(phased_pct)
        stats[f'{prefix}_h1'] = int(h1)
        stats[f'{prefix}_h1_pct'] = float(h1_pct)
        stats[f'{prefix}_h2'] = int(h2)
        stats[f'{prefix}_h2_pct'] = float(h2_pct)
        stats[f'{prefix}_htrans'] = int(htrans)
        stats[f'{prefix}_htrans_pct'] = float(htrans_pct)
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Parse and combine dip3d statistics files'
    )
    parser.add_argument(
        '--imputed',
        nargs='+',
        required=True,
        help='List of imputed stats files'
    )
    parser.add_argument(
        '--snp-tagged',
        nargs='+',
        required=True,
        help='List of snp-tagged stats files'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file path'
    )
    
    args = parser.parse_args()
    
    # Parse all stats files
    all_stats = []
    
    # Process imputed files
    for imputed_file in args.imputed:
        sample = Path(imputed_file).parts[-3]  # Extract sample name from path
        stats = parse_dip3d_stats(imputed_file, sample, 'imputed')
        all_stats.append(stats)
    
    # Process snp-tagged files
    for snp_file in args.snp_tagged:
        sample = Path(snp_file).parts[-3]  # Extract sample name from path
        stats = parse_dip3d_stats(snp_file, sample, 'snp-tagged')
        all_stats.append(stats)
    
    # Create DataFrame
    df = pd.DataFrame(all_stats)
    
    # Sort columns for better readability
    col_order = [
        'sample', 'tag_type', 
        'total_reads', 'reads_with_phased_frags', 'reads_with_phased_frags_pct',
        'total_frags', 'phased_frags', 'phased_frags_pct', 
        'h1_frags', 'h1_frags_pct', 'h2_frags', 'h2_frags_pct'
    ]
    
    # Add remaining columns (contact stats) in sorted order
    remaining_cols = [col for col in df.columns if col not in col_order]
    col_order.extend(sorted(remaining_cols))
    
    df = df[col_order]
    
    # Save to TSV
    df.to_csv(args.output, sep='\t', index=False)
    
    # Print summary
    print(f"Parsed {len(all_stats)} stats files")
    print(f"Samples: {', '.join(sorted(df['sample'].unique()))}")
    print(f"Output saved to: {args.output}")


if __name__ == '__main__':
    main()
