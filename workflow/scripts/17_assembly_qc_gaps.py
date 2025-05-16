#!/usr/bin/env python3
"""
Analyze N nucleotides in a FASTA file.
Produces a BED file with N locations and a statistics file.
"""

import argparse
import os
import re
from Bio import SeqIO
import statistics


def parse_args():
    parser = argparse.ArgumentParser(description='Analyze N nucleotides in a FASTA file')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output_prefix', required=True, help='Output prefix for files')
    return parser.parse_args()


def find_n_regions(sequence, seq_id):
    """Find regions of N nucleotides in a sequence"""
    n_regions = []
    n_pattern = re.compile(r'N+', re.IGNORECASE)
    
    for match in n_pattern.finditer(str(sequence)):
        start = match.start()
        end = match.end()
        n_regions.append((seq_id, start, end, end - start))
    
    return n_regions


def write_bed_file(n_regions, output_file):
    """Write N regions to a BED file"""
    with open(output_file, 'w') as f:
        for seq_id, start, end, length in n_regions:
            # BED format is 0-based for start and 1-based for end
            f.write(f"{seq_id}\t{start}\t{end}\tN_region\t{length}\n")


def calculate_stats(n_regions):
    """Calculate statistics for N regions"""
    if not n_regions:
        return {
            'total_n_count': 0,
            'total_gaps': 0,
            'avg_gap_size': 0,
            'gap_n50': 0,
            'max_gap_size': 0,
            'min_gap_size': 0
        }
    
    # Extract lengths of all N regions
    gap_lengths = [length for _, _, _, length in n_regions]
    total_n_count = sum(gap_lengths)
    total_gaps = len(gap_lengths)
    avg_gap_size = total_n_count / total_gaps if total_gaps > 0 else 0
    max_gap_size = max(gap_lengths) if gap_lengths else 0
    min_gap_size = min(gap_lengths) if gap_lengths else 0
    
    # Calculate N50
    gap_n50 = 0
    if gap_lengths:
        sorted_lengths = sorted(gap_lengths, reverse=True)
        cumulative_length = 0
        half_total = total_n_count / 2
        
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= half_total:
                gap_n50 = length
                break
    
    return {
        'total_n_count': total_n_count,
        'total_gaps': total_gaps,
        'avg_gap_size': avg_gap_size,
        'gap_n50': gap_n50,
        'max_gap_size': max_gap_size,
        'min_gap_size': min_gap_size
    }


def write_stats_file(stats, output_file):
    """Write statistics to a file"""
    with open(output_file, 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"Total N count\t{stats['total_n_count']}\n")
        f.write(f"Total gaps (N regions)\t{stats['total_gaps']}\n")
        f.write(f"Average gap size\t{stats['avg_gap_size']:.2f}\n")
        f.write(f"Gap N50\t{stats['gap_n50']}\n")
        f.write(f"Maximum gap size\t{stats['max_gap_size']}\n")
        f.write(f"Minimum gap size\t{stats['min_gap_size']}\n")


def main():
    args = parse_args()
    
    # Output file paths
    bed_file = f"{args.output_prefix}.n_regions.bed"
    stats_file = f"{args.output_prefix}.n_stats.tsv"
    
    all_n_regions = []
    sequence_stats = {}
    
    # Process each sequence in the FASTA file
    for record in SeqIO.parse(args.input, "fasta"):
        seq_id = record.id
        sequence = record.seq
        
        # Find N regions in this sequence
        n_regions = find_n_regions(sequence, seq_id)
        all_n_regions.extend(n_regions)
        
        # Calculate per-sequence statistics
        seq_stats = calculate_stats(n_regions)
        sequence_stats[seq_id] = seq_stats
    
    # Write BED file with all N regions
    write_bed_file(all_n_regions, bed_file)
    
    # Calculate and write overall statistics
    overall_stats = calculate_stats(all_n_regions)
    write_stats_file(overall_stats, stats_file)
    
    # Print summary to console
    print(f"Analysis complete for {args.input}")
    print(f"Found {overall_stats['total_gaps']} N regions with a total of {overall_stats['total_n_count']} N nucleotides")
    print(f"BED file written to: {bed_file}")
    print(f"Statistics written to: {stats_file}")


if __name__ == "__main__":
    main() 