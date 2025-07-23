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
    
    # Read PAF file - only read the first 12 columns we need
    colnames = ['read_name', 'read_length', 'read_start', 'read_end', 'strand', 
                'chrom', 'chrom_length', 'start', 'end', 'num_matches', 'alignment_length', 'mapq']
    
    try:
        paf_df = pd.read_csv(paf_file, header=None, sep="\t", names=colnames, usecols=range(12), low_memory=False)
    except Exception as e:
        print(f"Error reading PAF file: {e}")
        sys.exit(1)
    
    print(f"Read {len(paf_df)} alignments from PAF file")
    print(f"Unique chromosomes: {sorted(paf_df['chrom'].unique())}")
    
    # Calculate position (midpoint of alignment)
    try:
        # Convert to numeric, errors='coerce' will turn invalid values to NaN
        paf_df["start"] = pd.to_numeric(paf_df["start"], errors='coerce')
        paf_df["end"] = pd.to_numeric(paf_df["end"], errors='coerce')
        
        # Check for any NaN values after conversion
        if paf_df["start"].isna().any() or paf_df["end"].isna().any():
            print("Error: Non-numeric values found in 'start' or 'end' columns after conversion.")
            print("Rows with NaN in start column:")
            print(paf_df[paf_df["start"].isna()][["start", "end"]].head())
            print("Rows with NaN in end column:")
            print(paf_df[paf_df["end"].isna()][["start", "end"]].head())
            sys.exit(1)
            
    except Exception as e:
        print("Error: Failed to convert 'start' or 'end' columns to numeric.")
        print("start column head:\n", paf_df["start"].head())
        print("end column head:\n", paf_df["end"].head())
        raise e

    paf_df["position"] = ((paf_df["start"] + paf_df["end"]) / 2).astype(int)
    
    # Extract the base read name (remove the segment suffix)
    # e86e2934-4d83-467f-b750-8ea18c010bea_0001034108:000:0000000012:0115855612 -> e86e2934-4d83-467f-b750-8ea18c010bea_0001034108
    paf_df["base_read_name"] = paf_df["read_name"].str.split(":", n=1).str[0]
    
    print(f"Unique base read names: {paf_df['base_read_name'].nunique()}")

    # Filter for standard chromosomes (but include all chromosomes for now to debug)
    select_chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    print(f"Before chromosome filtering: {len(paf_df)} alignments")
    paf_df_filtered = paf_df[paf_df["chrom"].isin(select_chrs)]
    print(f"After chromosome filtering: {len(paf_df_filtered)} alignments")
    
    # If no alignments after filtering, use all chromosomes
    if len(paf_df_filtered) == 0:
        print("No alignments found for standard chromosomes, using all chromosomes")
        paf_df_filtered = paf_df
    
    # Group by base read name to find contacts within reads
    read_groups = paf_df_filtered.groupby("base_read_name")
    
    print(f"Found {len(read_groups)} unique base reads")
    
    # Check how many reads have multiple alignments
    multi_alignment_reads = sum(1 for name, group in read_groups if len(group) >= 2)
    print(f"Reads with multiple alignments: {multi_alignment_reads}")
    
    # Open output file and write header
    with open(pairs_file, 'w') as out_file:
        # Write pairs format header
        out_file.write("## pairs format v1.0\n")
        out_file.write("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 mapq1 mapq2\n")
        
        total_pairs = 0
        processed_reads = 0
        
        for base_read_name, group in read_groups:
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
                    out_file.write(f"{base_read_name}\t{row1['chrom']}\t{row1['position']}\t"
                                 f"{row2['chrom']}\t{row2['position']}\t{row1['strand']}\t"
                                 f"{row2['strand']}\t{row1['mapq']}\t{row2['mapq']}\n")
                    total_pairs += 1
            
            if processed_reads % 10000 == 0:
                print(f"Processed {processed_reads} reads, generated {total_pairs} pairs")
    
    print(f"Conversion complete. Processed {processed_reads} reads, generated {total_pairs} pairs")

if __name__ == "__main__":
    main()