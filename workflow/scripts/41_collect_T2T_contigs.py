#!/usr/bin/env python3
"""
Collect T2T contig information from motif, alignment, and gap files.
Match contigs between files and extract chromosome assignments plus gap statistics.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd

def read_motif_file(filepath):
    """Read motif file and return dictionary of contig -> length"""
    motif_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    contig = parts[0]
                    length = int(parts[1])
                    motif_data[contig] = length
    return motif_data

def read_alignment_file(filepath):
    """Read alignment file and return dictionary of contig -> chromosome info"""
    alignment_data = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    contig = parts[0]
                    contig_length = int(parts[1])
                    gap_length = int(parts[2])
                    chromosome = parts[3]
                    ref_length = int(parts[4])
                    
                    alignment_data[contig] = {
                        'contig_length': contig_length,
                        'gap_length': gap_length,
                        'chromosome': chromosome,
                        'ref_length': ref_length
                    }
    return alignment_data

def read_gap_bed_file(filepath):
    """Read gap BED file and return dictionary of contig -> gap statistics"""
    gap_data = {}
    
    if not filepath or not Path(filepath).exists():
        return gap_data
        
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    contig = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    gap_type = parts[3]  # Should be "N_region"
                    gap_size = int(parts[4])
                    
                    if contig not in gap_data:
                        gap_data[contig] = {
                            'n_gaps': 0,
                            'total_gap_length': 0,
                            'gap_sizes': [],
                            'gap_positions': []
                        }
                    
                    gap_data[contig]['n_gaps'] += 1
                    gap_data[contig]['total_gap_length'] += gap_size
                    gap_data[contig]['gap_sizes'].append(gap_size)
                    gap_data[contig]['gap_positions'].append((start, end))
    
    return gap_data

def collect_sample_data(sample, motif_file, alignment_file, gap_file):
    """Collect and match data for a single sample"""
    results = []
    
    if not Path(motif_file).exists():
        print(f"Warning: Motif file not found for {sample}: {motif_file}", file=sys.stderr)
        return results
    
    if not Path(alignment_file).exists():
        print(f"Warning: Alignment file not found for {sample}: {alignment_file}", file=sys.stderr)
        return results
    
    print(f"Processing sample: {sample}", file=sys.stderr)
    print(f"  Motif file: {motif_file}", file=sys.stderr)
    print(f"  Alignment file: {alignment_file}", file=sys.stderr)
    print(f"  Gap file: {gap_file}", file=sys.stderr)
    
    # Read all files
    motif_data = read_motif_file(motif_file)
    alignment_data = read_alignment_file(alignment_file)
    gap_data = read_gap_bed_file(gap_file)
    
    print(f"  Motif contigs: {len(motif_data)}", file=sys.stderr)
    print(f"  Alignment contigs: {len(alignment_data)}", file=sys.stderr)
    print(f"  Gap file contigs: {len(gap_data)}", file=sys.stderr)
    
    # Find contigs that appear in both motif and alignment files
    common_contigs = set(motif_data.keys()) & set(alignment_data.keys())
    print(f"  Common contigs (motif + alignment): {len(common_contigs)}", file=sys.stderr)
    
    # Process common contigs
    for contig in sorted(common_contigs):
        motif_length = motif_data[contig]
        alignment_info = alignment_data[contig]
        
        # Get gap information if available
        gap_info = gap_data.get(contig, {
            'n_gaps': 0,
            'total_gap_length': 0,
            'gap_sizes': [],
            'gap_positions': []
        })
        
        # Verify lengths match (with some tolerance for rounding)
        length_diff = abs(motif_length - alignment_info['contig_length'])
        if length_diff > 1000:  # Allow 1kb difference for rounding
            print(f"  Warning: Length mismatch for {contig}: motif={motif_length}, alignment={alignment_info['contig_length']}", file=sys.stderr)
        
        # Verify gap lengths match
        if gap_info['total_gap_length'] != alignment_info['gap_length']:
            print(f"  Warning: Gap length mismatch for {contig}: BED={gap_info['total_gap_length']}, alignment={alignment_info['gap_length']}", file=sys.stderr)
        
        # Calculate gap statistics
        mean_gap_size = gap_info['total_gap_length'] / gap_info['n_gaps'] if gap_info['n_gaps'] > 0 else 0
        max_gap_size = max(gap_info['gap_sizes']) if gap_info['gap_sizes'] else 0
        min_gap_size = min(gap_info['gap_sizes']) if gap_info['gap_sizes'] else 0
        
        results.append({
            'sample': sample,
            'contig': contig,
            'chromosome': alignment_info['chromosome'],
            'contig_length': alignment_info['contig_length'],
            'gap_length': alignment_info['gap_length'],
            'ref_length': alignment_info['ref_length'],
            'motif_length': motif_length,
            'n_gaps': gap_info['n_gaps'],
            'total_gap_length_bed': gap_info['total_gap_length'],
            'mean_gap_size': round(mean_gap_size, 1),
            'max_gap_size': max_gap_size,
            'min_gap_size': min_gap_size,
            'gap_density': round(gap_info['n_gaps'] / (alignment_info['contig_length'] / 1000000), 2) if alignment_info['contig_length'] > 0 else 0,  # gaps per Mb
            'classification': 'contig' if gap_info['n_gaps'] == 0 else 'scaffold'
        })
    
    # Report contigs only in specific files
    only_motif = set(motif_data.keys()) - set(alignment_data.keys())
    only_alignment = set(alignment_data.keys()) - set(motif_data.keys())
    only_gaps = set(gap_data.keys()) - common_contigs
    
    if only_motif:
        print(f"  Contigs only in motif file: {len(only_motif)}", file=sys.stderr)
    if only_alignment:
        print(f"  Contigs only in alignment file: {len(only_alignment)}", file=sys.stderr)
    if only_gaps:
        print(f"  Contigs only in gap file: {len(only_gaps)}", file=sys.stderr)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Collect T2T contig information including gap statistics')
    parser.add_argument('--motif-files', nargs='+', required=True, help='List of motif files')
    parser.add_argument('--alignment-files', nargs='+', required=True, help='List of alignment files')
    parser.add_argument('--gap-files', nargs='+', required=False, help='List of gap BED files (optional)')
    parser.add_argument('--samples', nargs='+', required=True, help='List of sample names (must match file order)')
    parser.add_argument('--output', required=True, help='Output file')
    
    args = parser.parse_args()
    
    # Validate input lengths
    if len(args.motif_files) != len(args.samples):
        print(f"Error: Number of motif files ({len(args.motif_files)}) must match number of samples ({len(args.samples)})", file=sys.stderr)
        sys.exit(1)
    
    if len(args.alignment_files) != len(args.samples):
        print(f"Error: Number of alignment files ({len(args.alignment_files)}) must match number of samples ({len(args.samples)})", file=sys.stderr)
        sys.exit(1)
    
    # Handle gap files - they might be optional or fewer than samples
    gap_files = args.gap_files if args.gap_files else [None] * len(args.samples)
    if args.gap_files and len(args.gap_files) != len(args.samples):
        print(f"Error: Number of gap files ({len(args.gap_files)}) must match number of samples ({len(args.samples)})", file=sys.stderr)
        sys.exit(1)
    
    all_results = []
    
    # Process each sample
    for i, sample in enumerate(args.samples):
        motif_file = args.motif_files[i]
        alignment_file = args.alignment_files[i]
        gap_file = gap_files[i] if gap_files[i] else None
        
        sample_results = collect_sample_data(sample, motif_file, alignment_file, gap_file)
        all_results.extend(sample_results)
    
    print(f"\nTotal collected records: {len(all_results)}", file=sys.stderr)
    
    # Create DataFrame and write output
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Sort by sample, then chromosome, then contig
        df = df.sort_values(['sample', 'chromosome', 'contig'])
        
        # Write to output file
        df.to_csv(args.output, sep='\t', index=False)
        
        # Print summary statistics
        print(f"\nSummary:", file=sys.stderr)
        print(f"  Total contigs: {len(df)}", file=sys.stderr)
        print(f"  Samples: {df['sample'].nunique()}", file=sys.stderr)
        print(f"  Chromosomes: {sorted(df['chromosome'].unique())}", file=sys.stderr)
        
        # Gap statistics
        total_gaps = df['n_gaps'].sum()
        contigs_with_gaps = (df['n_gaps'] > 0).sum()
        complete_contigs = (df['n_gaps'] == 0).sum()
        
        print(f"\nGap Statistics:", file=sys.stderr)
        print(f"  Total gaps across all contigs: {total_gaps}", file=sys.stderr)
        print(f"  Contigs with gaps: {contigs_with_gaps}", file=sys.stderr)
        print(f"  Complete contigs (no gaps): {complete_contigs}", file=sys.stderr)
        print(f"  Completion rate: {complete_contigs/len(df)*100:.1f}%", file=sys.stderr)
        
        # Per-sample statistics
        sample_stats = df.groupby('sample').agg({
            'contig': 'count',
            'n_gaps': 'sum',
            'gap_length': 'sum'
        }).rename(columns={'contig': 'n_contigs', 'n_gaps': 'total_gaps', 'gap_length': 'total_gap_bp'})
        
        print(f"\nPer-sample statistics:", file=sys.stderr)
        for sample, stats in sample_stats.iterrows():
            print(f"  {sample}: {stats['n_contigs']} contigs, {stats['total_gaps']} gaps, {stats['total_gap_bp']:,} bp gaps", file=sys.stderr)
            
    else:
        print("No data collected - creating empty output file", file=sys.stderr)
        # Create empty file with headers
        with open(args.output, 'w') as f:
            f.write("sample\tcontig\tchromosome\tcontig_length\tgap_length\tref_length\tmotif_length\tn_gaps\ttotal_gap_length_bed\tmean_gap_size\tmax_gap_size\tmin_gap_size\tgap_density\tclassification\n")

if __name__ == '__main__':
    main()