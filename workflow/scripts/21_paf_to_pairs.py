#!/usr/bin/env python
# coding: utf-8

# Convert cleaned PAF file to pairs format for Hi-C analysis
import sys
import pandas as pd
import itertools

def main():
    if len(sys.argv) != 3:
        print("Usage: python paf_to_pairs.py <input.clean.paf> <output.pairs>")
        sys.exit(1)
    
    paf_file = sys.argv[1]
    pairs_file = sys.argv[2]
    
    print(f"Converting PAF file {paf_file} to pairs format")
    
    # Read PAF file
    colnames = ['read_name', 'read_length', 'read_start', 'read_end', 'strand', 
                'chrom', 'chrom_length', 'start', 'end', 'num_matches', 'alignment_length', 'mapq']
    
    try:
        paf_df = pd.read_csv(paf_file, header=None, sep="\t", names=colnames, low_memory=False)
    except Exception as e:
        print(f"Error reading PAF file: {e}")
        sys.exit(1)
    
    # Calculate position (midpoint of alignment)
    paf_df["position"] = ((paf_df["start"] + paf_df["end"]) / 2).astype(int)
    
    # Filter for standard chromosomes
    select_chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    paf_df = paf_df[paf_df["chrom"].isin(select_chrs)]
    
    # Group by read name to find contacts within reads
    read_groups = paf_df.groupby("read_name")
    
    # Open output file and write header
    with open(pairs_file, 'w') as out_file:
        # Write pairs format header
        out_file.write("## pairs format v1.0\n")
        out_file.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 mapq1 mapq2\n")
        
        total_pairs = 0
        processed_reads = 0
        
        for read_name, group in read_groups:
            processed_reads += 1
            
            if len(group) < 2:
                continue
                
            # Sort by chromosome and position
            group = group.sort_values(['chrom', 'position']).reset_index(drop=True)
            
            # Generate all pairwise combinations
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    row1 = group.iloc[i]
                    row2 = group.iloc[j]
                    
                    # Write pair
                    out_file.write(f"{read_name}\t{row1['chrom']}\t{row1['position']}\t"
                                 f"{row2['chrom']}\t{row2['position']}\t{row1['strand']}\t"
                                 f"{row2['strand']}\t{row1['mapq']}\t{row2['mapq']}\n")
                    total_pairs += 1
            
            if processed_reads % 10000 == 0:
                print(f"Processed {processed_reads} reads, generated {total_pairs} pairs")
    
    print(f"Conversion complete. Processed {processed_reads} reads, generated {total_pairs} pairs")

if __name__ == "__main__":
    main()