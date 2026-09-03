#!/usr/bin/env python3
"""
Parse yak trioeval phase blocks (F records) for a single assembly and annotate with reference chromosomes.

This script processes both haplotypes of an assembly:
1. Extracts phase blocks from trioeval output files
2. Maps them to reference chromosomes using PAF alignments
3. Calculates switch positions
4. Outputs annotated phase blocks and switches for the assembly
"""

import pandas as pd
import argparse
from pathlib import Path

def parse_paf_file(paf_path):
    """
    Parse PAF file to create contig-to-chromosome mapping.
    Returns dict: contig_name -> (chr, ref_start, ref_end, strand, mapq, query_start, query_end)
    """
    mappings = {}
    
    with open(paf_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            query_name = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            num_matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])
            
            # Only keep primary alignments with good quality
            if mapq >= 5:
                # Store best mapping for each contig (highest MAPQ)
                if query_name not in mappings or mappings[query_name][4] < mapq:
                    mappings[query_name] = (
                        target_name,
                        target_start,
                        target_end,
                        strand,
                        mapq,
                        query_start,
                        query_end
                    )
    
    return mappings


def strip_pansn(name):
    """Strip PanSN prefix ({sample}#{hap}#{contig}) returning the bare contig name."""
    parts = name.split('#')
    return parts[2] if len(parts) == 3 else name


def parse_trioeval_phaseblocks(filepath, contig_mappings):
    """Parse phase blocks (F records) from trioeval file and annotate with ref positions."""

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
    
    phaseblocks = []
    
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            if not fields or fields[0] != 'F':
                continue
            
            contig_name = fields[1]
            block_type = int(fields[2])  # 1=paternal, 2=maternal
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            kmer_count = int(fields[5])

            # PanSN-named FASTAs (pofo assemblies) carry sample#hap#contig headers;
            # strip prefix so the bare contig name matches the PAF file.
            lookup_name = strip_pansn(contig_name)

            # Get chromosome mapping
            if lookup_name in contig_mappings:
                chr_name, ref_start, ref_end, strand, mapq, query_start, query_end = contig_mappings[lookup_name]
                
                # Calculate reference positions
                # Simple linear interpolation
                contig_len = query_end - query_start
                if contig_len > 0:
                    ratio_start = (start_pos - query_start) / contig_len
                    ratio_end = (end_pos - query_start) / contig_len
                    
                    if strand == '+':
                        block_ref_start = int(ref_start + ratio_start * (ref_end - ref_start))
                        block_ref_end = int(ref_start + ratio_end * (ref_end - ref_start))
                    else:
                        # Reverse strand
                        block_ref_start = int(ref_end - ratio_end * (ref_end - ref_start))
                        block_ref_end = int(ref_end - ratio_start * (ref_end - ref_start))
                else:
                    block_ref_start = ref_start
                    block_ref_end = ref_end
            else:
                chr_name = "unmapped"
                block_ref_start = None
                block_ref_end = None
                strand = None
                mapq = 0
            
            phaseblocks.append({
                'sample': sample_base,
                'method': method,
                'haplotype': haplotype,
                'contig': contig_name,
                'block_type': 'paternal' if block_type == 1 else 'maternal',
                'start': start_pos,
                'end': end_pos,
                'length': end_pos - start_pos,
                'kmer_count': kmer_count,
                'chr': chr_name,
                'ref_start': block_ref_start,
                'ref_end': block_ref_end,
                'strand': strand,
                'mapq': mapq
            })
    
    return pd.DataFrame(phaseblocks)


def calculate_switch_positions(phaseblocks_df, min_kmer_support=5):
    """
    Identify switch positions from phase blocks.
    
    Args:
        phaseblocks_df: DataFrame with phase blocks
        min_kmer_support: Minimum k-mer count for both blocks at a transition
    """
    
    switches = []
    
    # Group by contig (across both haplotypes)
    for (sample, method, haplotype, contig), group in phaseblocks_df.groupby(
        ['sample', 'method', 'haplotype', 'contig']
    ):
        contig_blocks = group.sort_values('start').reset_index(drop=True)
        
        for i in range(len(contig_blocks) - 1):
            curr_block = contig_blocks.iloc[i]
            next_block = contig_blocks.iloc[i + 1]
            
            # Check if types differ (= switch)
            if curr_block['block_type'] != next_block['block_type']:
                # Filter out low-confidence switches
                if (curr_block['kmer_count'] >= min_kmer_support and 
                    next_block['kmer_count'] >= min_kmer_support):
                    
                    # Calculate reference position of switch
                    switch_ref_pos = None
                    if pd.notna(next_block['ref_start']):
                        switch_ref_pos = next_block['ref_start']
                    
                    switches.append({
                        'sample': sample,
                        'method': method,
                        'haplotype': haplotype,
                        'contig': contig,
                        'switch_position': next_block['start'],
                        'from_type': curr_block['block_type'],
                        'to_type': next_block['block_type'],
                        'from_kmer_count': curr_block['kmer_count'],
                        'to_kmer_count': next_block['kmer_count'],
                        'from_block_length': curr_block['length'],
                        'to_block_length': next_block['length'],
                        'chr': next_block['chr'],
                        'ref_position': switch_ref_pos,
                        'strand': next_block['strand']
                    })
    
    return pd.DataFrame(switches)


def main():
    """Main function."""
    
    parser = argparse.ArgumentParser(
        description='Parse yak trioeval phase blocks for a single assembly and annotate with reference positions.'
    )
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        help='Input trioeval files (both haplotypes)'
    )
    parser.add_argument(
        '--paf',
        nargs='+',
        required=True,
        help='PAF alignment files (both haplotypes)'
    )
    parser.add_argument(
        '-p', '--phaseblocks',
        required=True,
        help='Output phase blocks TSV file'
    )
    parser.add_argument(
        '-w', '--switches',
        required=True,
        help='Output switch positions TSV file'
    )
    parser.add_argument(
        '--min-kmer-support',
        type=int,
        default=5,
        help='Minimum k-mer support for both blocks at a switch (default: 5)'
    )
    
    args = parser.parse_args()
    
    if len(args.input) != 2:
        raise ValueError(f"Expected 2 trioeval files (haplotype1 and haplotype2), got {len(args.input)}")
    
    if len(args.paf) != 2:
        raise ValueError(f"Expected 2 PAF files (haplotype1 and haplotype2), got {len(args.paf)}")
    
    min_kmer_support = args.min_kmer_support
    
    print(f"Processing assembly with {len(args.input)} haplotypes...")
    print(f"Minimum k-mer support for switches: {min_kmer_support}")
    
    # Load PAF mappings for each haplotype
    all_mappings = {}
    for paf_path in args.paf:
        paf_path = Path(paf_path)
        print(f"  Loading PAF: {paf_path}")
        
        # Extract haplotype from filename
        if 'haplotype1' in str(paf_path):
            hap = 'haplotype1'
        elif 'haplotype2' in str(paf_path):
            hap = 'haplotype2'
        else:
            print(f"  Warning: Could not determine haplotype from {paf_path}")
            continue
        
        mappings = parse_paf_file(paf_path)
        all_mappings[hap] = mappings
        print(f"    Loaded {len(mappings)} contig mappings for {hap}")

    # Merge all PAF mappings into one lookup so that contigs reassigned across
    # haplotypes by pofo phasing (hap2 contigs appearing in hap1 trioeval) are
    # still resolved. Contig names are unique across haplotypes, so merging is safe.
    merged_mappings = {}
    for m in all_mappings.values():
        merged_mappings.update(m)

    # Parse phase blocks for both haplotypes
    all_phaseblocks = []

    for filepath in args.input:
        filepath = Path(filepath)
        print(f"  Processing: {filepath}")

        # Determine haplotype
        if 'haplotype1' in str(filepath):
            hap = 'haplotype1'
        elif 'haplotype2' in str(filepath):
            hap = 'haplotype2'
        else:
            print(f"  Warning: Could not determine haplotype from {filepath}")
            continue

        contig_mappings = merged_mappings
        
        phaseblocks_df = parse_trioeval_phaseblocks(filepath, contig_mappings)
        
        if not phaseblocks_df.empty:
            all_phaseblocks.append(phaseblocks_df)
            print(f"    Found {len(phaseblocks_df)} phase blocks")
            mapped = (phaseblocks_df['chr'] != 'unmapped').sum()
            print(f"    Mapped to reference: {mapped}/{len(phaseblocks_df)} ({100*mapped/len(phaseblocks_df):.1f}%)")
    
    # Combine and save phase blocks
    if all_phaseblocks:
        phaseblocks_df = pd.concat(all_phaseblocks, ignore_index=True)
        phaseblocks_df.to_csv(args.phaseblocks, sep='\t', index=False)
        print(f"\nSaved phase blocks data: {len(phaseblocks_df)} blocks")
        
        # Print per-haplotype statistics
        print("\n=== Phase Block Statistics ===")
        pb_stats = phaseblocks_df.groupby(['haplotype', 'block_type']).agg({
            'length': ['count', 'sum', 'mean'],
            'kmer_count': ['mean', 'median']
        }).round(2)
        print(pb_stats)
        
        # Calculate and save switches
        switches_df = calculate_switch_positions(phaseblocks_df, min_kmer_support)
        switches_df.to_csv(args.switches, sep='\t', index=False)
        print(f"\nSaved switch positions: {len(switches_df)} switches (filtered with min k-mer support={min_kmer_support})")
        
        # Print statistics
        if not switches_df.empty:
            print("\n=== Switch Statistics ===")
            switch_stats = switches_df.groupby('haplotype').agg({
                'switch_position': 'count',
                'from_kmer_count': ['mean', 'median'],
                'to_kmer_count': ['mean', 'median']
            }).round(2)
            switch_stats.columns = ['n_switches', 'mean_from_kmers', 'median_from_kmers', 
                                   'mean_to_kmers', 'median_to_kmers']
            print(switch_stats)
            
            # Print per-chromosome statistics
            print("\n=== Switches per Chromosome ===")
            chr_stats = switches_df[switches_df['chr'] != 'unmapped'].groupby('chr').size().sort_values(ascending=False)
            print(chr_stats)
    else:
        print("\nWarning: No phase blocks found!")
        pd.DataFrame().to_csv(args.phaseblocks, sep='\t', index=False)
        pd.DataFrame().to_csv(args.switches, sep='\t', index=False)

if __name__ == '__main__':
    main()