#!/usr/bin/env python3
"""
Merge multiple BED files from multiple samples into a single aggregated table.
Extracts haplotype and contig information from contig names.
"""

import argparse
import sys
import pandas as pd
import re

def parse_contig_name(contig):
    """Parse contig name to extract haplotype and contig ID."""
    # Expected format: haplotype1-0000001, haplotype2-0000001, unassigned-0000001
    match = re.match(r'^(haplotype[12]|unassigned)-(.+)$', contig)
    if match:
        haplotype = match.group(1)
        chr_id = match.group(2)
        return haplotype, chr_id
    else:
        # If pattern doesn't match, return original as chr and unknown haplotype
        return 'unknown', contig

def read_bed_file(bed_file, tool_name, sample_name):
    """Read BED file and return as list of dictionaries."""
    regions = []
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 3:
                    continue
                
                # Parse haplotype and chr from contig name
                contig = fields[0]
                haplotype, chr_id = parse_contig_name(contig)
                
                # Parse standard BED fields
                region = {
                    'sample': sample_name,
                    'source': tool_name,
                    'haplotype': haplotype,
                    'chr': chr_id,
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'feature_type': fields[3] if len(fields) > 3 else '.',
                }
                
                regions.append(region)
    except Exception as e:
        print(f"WARNING: Error reading {bed_file}: {e}", file=sys.stderr)
    
    return regions

def merge_bed_files(flagger_files, nucflag_files, gaps_files, samples, output_file):
    """Merge BED files from different tools and samples into single table."""
    
    all_regions = []
    
    for i, sample in enumerate(samples):
        print(f"\nProcessing sample: {sample}", file=sys.stderr)
        
        # Read Flagger
        print(f"  Reading Flagger: {flagger_files[i]}", file=sys.stderr)
        flagger_regions = read_bed_file(flagger_files[i], 'Flagger', sample)
        print(f"    Found {len(flagger_regions)} regions", file=sys.stderr)
        all_regions.extend(flagger_regions)
        
        # Read Nucflag
        print(f"  Reading Nucflag: {nucflag_files[i]}", file=sys.stderr)
        nucflag_regions = read_bed_file(nucflag_files[i], 'Nucflag', sample)
        print(f"    Found {len(nucflag_regions)} regions", file=sys.stderr)
        all_regions.extend(nucflag_regions)
        
        # Read Gaps
        print(f"  Reading Gaps: {gaps_files[i]}", file=sys.stderr)
        gaps_regions = read_bed_file(gaps_files[i], 'Gaps', sample)
        print(f"    Found {len(gaps_regions)} regions", file=sys.stderr)
        all_regions.extend(gaps_regions)
    
    if not all_regions:
        print("ERROR: No regions found!", file=sys.stderr)
        sys.exit(1)
    
    # Convert to DataFrame
    df = pd.DataFrame(all_regions)
    
    # Ensure column order
    columns = ['sample', 'source', 'haplotype', 'chr', 'start', 'end', 'feature_type']
    df = df[columns]
    
    # Sort by sample, haplotype, chromosome, and position
    df = df.sort_values(['sample', 'haplotype', 'chr', 'start', 'end'])
    
    # Write output
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\nSUMMARY:", file=sys.stderr)
    print(f"Total samples: {df['sample'].nunique()}", file=sys.stderr)
    print(f"Total regions: {len(df)}", file=sys.stderr)
    print(f"\nRegions by source:", file=sys.stderr)
    for source, count in df['source'].value_counts().items():
        print(f"  {source}: {count}", file=sys.stderr)
    print(f"\nRegions by haplotype:", file=sys.stderr)
    for haplotype, count in df['haplotype'].value_counts().items():
        print(f"  {haplotype}: {count}", file=sys.stderr)
    print(f"\nRegions by sample:", file=sys.stderr)
    for sample, count in df.groupby('sample').size().items():
        print(f"  {sample}: {count}", file=sys.stderr)
    print(f"\nOutput written to: {output_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Merge BED files from multiple samples into single table')
    parser.add_argument('--flagger', nargs='+', required=True, help='Flagger BED files')
    parser.add_argument('--nucflag', nargs='+', required=True, help='Nucflag BED files')
    parser.add_argument('--gaps', nargs='+', required=True, help='Gaps BED files')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample names')
    parser.add_argument('--output', required=True, help='Output TSV file')
    
    args = parser.parse_args()
    
    # Verify all lists have same length
    if not (len(args.flagger) == len(args.nucflag) == len(args.gaps) == len(args.samples)):
        print("ERROR: Number of flagger, nucflag, gaps files and samples must match!", file=sys.stderr)
        sys.exit(1)
    
    merge_bed_files(args.flagger, args.nucflag, args.gaps, args.samples, args.output)

if __name__ == '__main__':
    main()