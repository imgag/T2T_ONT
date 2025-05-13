#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
import uuid


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess diploid assembly for rfhap')
    parser.add_argument('-a', '--hap1', required=True, help='Path to haplotype 1 FASTA file')
    parser.add_argument('-b', '--hap2', required=True, help='Path to haplotype 2 FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for files')
    parser.add_argument('-c', '--chunk_size', type=int, default=10000, help='Chunk size in bp (default: 10000)')
    parser.add_argument('-p', '--prefix1', default='hapA_', help='Prefix for haplotype 1 chunks')
    parser.add_argument('-q', '--prefix2', default='hapB_', help='Prefix for haplotype 2 chunks')
    return parser.parse_args()


def split_fasta_into_chunks(fasta_file, chunk_size, prefix):
    """Split a FASTA file into chunks of specified size"""
    chunks = []
    chunk_positions = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        contig_id = record.id
        contig_len = len(seq)
        
        for i in range(0, contig_len, chunk_size):
            end = min(i + chunk_size, contig_len)
            chunk_seq = seq[i:end]
            
            # Create unique chunk ID that encodes original position
            chunk_id = f"{prefix}{contig_id}:{i+1}-{end}"
            
            # Store chunk info for GFA creation
            chunks.append((chunk_id, chunk_seq))
            chunk_positions.append((chunk_id, contig_id, i+1, end))
    
    return chunks, chunk_positions


def write_merged_fasta(chunks_hap1, chunks_hap2, output_file):
    """Write all chunks to a single FASTA file"""
    with open(output_file, 'w') as f:
        for chunk_id, chunk_seq in chunks_hap1 + chunks_hap2:
            f.write(f">{chunk_id}\n")
            f.write(f"{chunk_seq}\n")


def create_gfa(chunks_positions_hap1, chunks_positions_hap2, output_file):
    """Create a GFA file where chunks from the same contig are linked"""
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write nodes (segments)
        all_chunks = chunks_positions_hap1 + chunks_positions_hap2
        for chunk_id, contig_id, start, end in all_chunks:
            length = end - start + 1
            f.write(f"S\t{chunk_id}\t*\tLN:i:{length}\tOR:Z:{contig_id}\tST:i:{start}\tEN:i:{end}\n")
        
        # Write edges (links)
        # Group chunks by contig
        contigs = {}
        for chunk_id, contig_id, start, end in all_chunks:
            if contig_id not in contigs:
                contigs[contig_id] = []
            contigs[contig_id].append((chunk_id, start, end))
        
        # Sort chunks by position and create links
        for contig_id, chunks in contigs.items():
            sorted_chunks = sorted(chunks, key=lambda x: x[1])
            for i in range(len(sorted_chunks) - 1):
                chunk1_id = sorted_chunks[i][0]
                chunk2_id = sorted_chunks[i + 1][0]
                f.write(f"L\t{chunk1_id}\t+\t{chunk2_id}\t+\t0M\n")


def create_chunk_info(chunks_positions_hap1, chunks_positions_hap2, output_file):
    """Create a TSV file with chunk information for traceability"""
    with open(output_file, 'w') as f:
        f.write("chunk_id\tcontig_id\tstart\tend\thaplotype\n")
        
        for chunk_id, contig_id, start, end in chunks_positions_hap1:
            f.write(f"{chunk_id}\t{contig_id}\t{start}\t{end}\thap1\n")
            
        for chunk_id, contig_id, start, end in chunks_positions_hap2:
            f.write(f"{chunk_id}\t{contig_id}\t{start}\t{end}\thap2\n")


def main():
    args = parse_args()
    
    # Split haplotype 1 into chunks
    chunks_hap1, chunks_positions_hap1 = split_fasta_into_chunks(
        args.hap1, args.chunk_size, args.prefix1
    )
    
    # Split haplotype 2 into chunks
    chunks_hap2, chunks_positions_hap2 = split_fasta_into_chunks(
        args.hap2, args.chunk_size, args.prefix2
    )
    
    # Write merged FASTA
    merged_fasta = f"{args.output}.fasta"
    write_merged_fasta(chunks_hap1, chunks_hap2, merged_fasta)
    
    # Create GFA
    gfa_file = f"{args.output}.gfa"
    create_gfa(chunks_positions_hap1, chunks_positions_hap2, gfa_file)
    
    # Create chunk info file for traceability
    chunk_info = f"{args.output}.chunks.tsv"
    create_chunk_info(chunks_positions_hap1, chunks_positions_hap2, chunk_info)
    
    print(f"Processed {len(chunks_hap1)} chunks from haplotype 1")
    print(f"Processed {len(chunks_hap2)} chunks from haplotype 2")
    print(f"Merged FASTA written to: {merged_fasta}")
    print(f"GFA written to: {gfa_file}")
    print(f"Chunk information written to: {chunk_info}")


if __name__ == "__main__":
    main()