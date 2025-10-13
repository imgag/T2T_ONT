#!/usr/bin/env python3
"""
Combine per-chromosome RFMix MSP files into a single genome-wide file.
"""

import pandas as pd
import argparse
import sys

def combine_msp_files(msp_files, output_file):
    """Combine multiple MSP files into one."""
    
    dfs = []
    header_lines = []
    
    for msp_file in msp_files:
        # Read the file to extract header lines
        with open(msp_file, 'r') as f:
            lines = f.readlines()
        
        # Find header lines
        subpop_line = None
        col_header_line = None
        data_start_idx = 0
        
        for i, line in enumerate(lines):
            if line.startswith('#Subpopulation'):
                if not header_lines:  # Keep only from first file
                    subpop_line = line
            elif line.startswith('#chm'):
                col_header_line = line
                data_start_idx = i + 1
                break
        
        # Store header lines from first file only
        if not header_lines and subpop_line:
            header_lines.append(subpop_line)
        
        # Parse column names from the header line
        if col_header_line:
            col_names = col_header_line.strip().split('\t')
            # Remove '#' from first column name
            col_names[0] = col_names[0].lstrip('#')
        else:
            raise ValueError(f"Could not find column header in {msp_file}")
        
        # Read data starting from the correct line
        df = pd.read_csv(msp_file, sep='\t', skiprows=data_start_idx, 
                        names=col_names, header=None)
        
        print(f"Read {msp_file}: {len(df)} segments, columns: {col_names[:5]}", 
              file=sys.stderr)
        
        dfs.append(df)
    
    # Combine all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)
    
    print(f"Combined dataframe: {len(combined_df)} segments", file=sys.stderr)
    print(f"Columns: {combined_df.columns.tolist()[:8]}", file=sys.stderr)
    
    # Get the chromosome column name (first column, should be 'chm')
    chrom_col = combined_df.columns[0]
    pos_col = combined_df.columns[1]  # 'spos'
    
    # Sort by chromosome and position
    # Extract numeric part of chromosome for sorting
    combined_df['chr_num'] = (combined_df[chrom_col]
                             .str.replace('chr', '', regex=False)
                             .replace('X', '23')
                             .replace('Y', '24')
                             .astype(int))
    
    combined_df = combined_df.sort_values(['chr_num', pos_col]).drop('chr_num', axis=1)
    
    print(f"Sorted by chromosome and position", file=sys.stderr)
    
    # Write output
    with open(output_file, 'w') as f:
        # Write header lines
        for line in header_lines:
            f.write(line)
        
        # Write column header (restore the # prefix for first column)
        cols = ['#' + combined_df.columns[0]] + list(combined_df.columns[1:])
        f.write('\t'.join(cols) + '\n')
        
        # Write data
        combined_df.to_csv(f, sep='\t', index=False, header=False)
    
    print(f"Combined MSP file written to {output_file} with {len(combined_df)} total segments", 
          file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Combine per-chromosome RFMix MSP files'
    )
    parser.add_argument('--msp-files', nargs='+', required=True,
                       help='Input RFMix MSP files (one per chromosome)')
    parser.add_argument('--output', required=True,
                       help='Output combined MSP file')
    
    args = parser.parse_args()
    
    combine_msp_files(args.msp_files, args.output)

if __name__ == '__main__':
    main()