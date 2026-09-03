#!/usr/bin/env python3
"""
Parse yak trioeval output files and create summary tables.

Output format (per sequence):
S  seqName     #patKmer  #matKmer  #pat-pat  #pat-mat  #mat-pat  #mat-mat  seqLen

Summary format (per assembly):
W  #switchErr  denominator  switchErrRate
H  #hammingErr denominator  hammingErrRate
N  #totPatKmer #totMatKmer  errRate
"""

import pandas as pd
import argparse
from pathlib import Path

def parse_trioeval_file(filepath, min_length=500000):
    """Parse a single yak trioeval output file for S, W, H, N records only."""
    
    sample = Path(filepath).stem.replace('.trio_eval', '')
    
    # Extract sample name and haplotype
    if '.haplotype1' in sample:
        haplotype = 'haplotype1'
        sample = sample.replace('.haplotype1', '')
    elif '.haplotype2' in sample:
        haplotype = 'haplotype2'
        sample = sample.replace('.haplotype2', '')
    elif '.maternal' in sample:
        haplotype = 'maternal'
        sample = sample.replace('.maternal', '')
    elif '.paternal' in sample:
        haplotype = 'paternal'
        sample = sample.replace('.paternal', '')
    else:
        haplotype = 'unknown'

    # Determine phasing method and clean sample name
    if 'hifiasm_trio' in sample:
        method = 'hifiasm_trio'
        sample_base = sample.replace('_hifiasm_trio', '')
    elif 'hifiasm_porec' in sample:
        method = 'hifiasm_porec'
        sample_base = sample.replace('_hifiasm_porec', '')
    elif '_trio' in sample:
        method = 'verkko_trio'
        sample_base = sample.replace('_trio', '')
    elif '_pofo' in sample:
        method = 'verkko_pofo'
        sample_base = sample.replace('_pofo', '')
    else:
        method = 'verkko_porec'
        sample_base = sample
    
    contigs = []
    summary = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            if not fields:
                continue
                
            record_type = fields[0]
            
            # Per-sequence data
            if record_type == 'S':
                contig_name = fields[1]
                pat_kmer = int(fields[2])
                mat_kmer = int(fields[3])
                pat_pat = int(fields[4])
                pat_mat = int(fields[5])
                mat_pat = int(fields[6])
                mat_mat = int(fields[7])
                seq_len = int(fields[8])
                
                # Filter by length
                if seq_len < min_length:
                    continue
                
                # Calculate per-contig metrics
                total_kmers = pat_kmer + mat_kmer
                if total_kmers > 0:
                    pat_fraction = pat_kmer / total_kmers
                    mat_fraction = mat_kmer / total_kmers
                else:
                    pat_fraction = 0
                    mat_fraction = 0
                
                # Calculate phasing consistency
                if pat_kmer > 0:
                    pat_consistency = pat_pat / pat_kmer
                else:
                    pat_consistency = 0
                    
                if mat_kmer > 0:
                    mat_consistency = mat_mat / mat_kmer
                else:
                    mat_consistency = 0
                
                # Determine likely haplotype origin
                if pat_fraction > 0.7:
                    predicted_haplotype = 'paternal'
                elif mat_fraction > 0.7:
                    predicted_haplotype = 'maternal'
                else:
                    predicted_haplotype = 'mixed'
                
                contigs.append({
                    'sample': sample_base,
                    'method': method,
                    'haplotype': haplotype,
                    'contig': contig_name,
                    'length': seq_len,
                    'pat_kmer': pat_kmer,
                    'mat_kmer': mat_kmer,
                    'pat_pat': pat_pat,
                    'pat_mat': pat_mat,
                    'mat_pat': mat_pat,
                    'mat_mat': mat_mat,
                    'total_kmers': total_kmers,
                    'pat_fraction': pat_fraction,
                    'mat_fraction': mat_fraction,
                    'pat_consistency': pat_consistency,
                    'mat_consistency': mat_consistency,
                    'predicted_origin': predicted_haplotype
                })
            
            # Assembly-wide metrics
            elif record_type == 'W':  # Switch errors
                summary['switch_errors'] = int(fields[1])
                summary['switch_denominator'] = int(fields[2])
                summary['switch_error_rate'] = float(fields[3])
                
            elif record_type == 'H':  # Hamming errors
                summary['hamming_errors'] = int(fields[1])
                summary['hamming_denominator'] = int(fields[2])
                summary['hamming_error_rate'] = float(fields[3])
                
            elif record_type == 'N':  # Total k-mers
                summary['total_pat_kmer'] = int(fields[1])
                summary['total_mat_kmer'] = int(fields[2])
                summary['total_error_rate'] = float(fields[3])
    
    # Convert to DataFrames
    contigs_df = pd.DataFrame(contigs) if contigs else pd.DataFrame()
    
    # Add summary info to each row
    if not contigs_df.empty:
        for key, value in summary.items():
            contigs_df[f'asm_{key}'] = value
    
    return contigs_df, summary, sample_base, method, haplotype


def main():
    """Main function."""
    
    parser = argparse.ArgumentParser(
        description='Parse yak trioeval output files and create summary tables.'
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input trioeval files'
    )
    parser.add_argument(
        '-s', '--summary',
        required=True,
        help='Output summary TSV file'
    )
    parser.add_argument(
        '-c', '--per-contig',
        required=True,
        help='Output per-contig TSV file'
    )
    parser.add_argument(
        '-m', '--min-length',
        type=int,
        default=500000,
        help='Minimum contig length (default: 500000)'
    )
    
    args = parser.parse_args()
    
    min_length = args.min_length
    
    all_contigs = []
    all_summaries = []
    
    print(f"Processing {len(args.input)} trioeval files...")
    print(f"Minimum contig length: {min_length:,} bp")
    
    for filepath in args.input:
        print(f"  Processing: {filepath}")
        
        contigs_df, summary, sample, method, haplotype = parse_trioeval_file(
            filepath, min_length
        )
        
        if not contigs_df.empty:
            all_contigs.append(contigs_df)
            print(f"    Found {len(contigs_df)} contigs >= {min_length:,} bp")
        
        # Add sample info to summary
        summary['sample'] = sample
        summary['method'] = method
        summary['haplotype'] = haplotype
        all_summaries.append(summary)
    
    # Combine all data
    if all_contigs:
        per_contig_df = pd.concat(all_contigs, ignore_index=True)
        per_contig_df.to_csv(args.per_contig, sep='\t', index=False)
        print(f"\nSaved per-contig data: {len(per_contig_df)} rows")
    else:
        print("\nWarning: No contigs found!")
        pd.DataFrame().to_csv(args.per_contig, sep='\t', index=False)
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(all_summaries)
    
    # Reorder columns
    col_order = ['sample', 'method', 'haplotype', 
                 'switch_error_rate', 'switch_errors', 'switch_denominator',
                 'hamming_error_rate', 'hamming_errors', 'hamming_denominator',
                 'total_error_rate', 'total_pat_kmer', 'total_mat_kmer']
    
    summary_df = summary_df[col_order]
    summary_df.to_csv(args.summary, sep='\t', index=False)
    print(f"Saved summary data: {len(summary_df)} assemblies")
    
    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(summary_df.to_string(index=False))
    
    if all_contigs:
        print("\n=== Per-Haplotype Contig Statistics ===")
        stats = per_contig_df.groupby(['sample', 'method', 'haplotype']).agg({
            'length': ['count', 'sum', 'mean', 'median'],
            'pat_fraction': 'mean',
            'mat_fraction': 'mean',
            'predicted_origin': lambda x: x.value_counts().to_dict()
        }).round(4)
        print(stats)
        
if __name__ == '__main__':
    main()