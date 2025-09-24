#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import glob

def main():
    parser = argparse.ArgumentParser(description='Create iAdmix frequency file from PLINK data')
    parser.add_argument('--plink-freq', required=True, nargs='+', help='PLINK frequency files (multiple files)')
    parser.add_argument('--pvar', required=True, help='PVAR file with variant info')
    parser.add_argument('--psam', required=True, help='PSAM file with population info')
    parser.add_argument('--output', required=True, help='Output frequency file for iAdmix')
    
    args = parser.parse_args()
    
    print(f"Processing {len(args.plink_freq)} frequency files:")
    for f in args.plink_freq:
        print(f"  - {f}")
    
    # Read PVAR file for variant information (handle comments properly)
    print("Reading PVAR file...")
    try:
        # PVAR files can have comments starting with ##
        # First, find the actual header line
        with open(args.pvar, 'r') as f:
            lines = f.readlines()
        
        # Find the header line (starts with #CHROM or CHROM)
        header_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('#CHROM') or (line.startswith('CHROM') and not line.startswith('##')):
                header_idx = i
                break
        
        print(f"Found header at line {header_idx + 1}")
        
        # Read the file starting from the header
        pvar_df = pd.read_csv(args.pvar, sep='\t', skiprows=header_idx, comment='#')
        
        # Clean column names
        pvar_df.columns = pvar_df.columns.str.replace('#', '')
        
    except Exception as e:
        print(f"Error reading PVAR file: {e}")
        print("Trying alternative parsing method...")
        
        # Alternative: read with comment='#' and infer header
        try:
            pvar_df = pd.read_csv(args.pvar, sep='\t', comment='#')
            if pvar_df.columns[0].startswith('#'):
                pvar_df.columns = pvar_df.columns.str.replace('#', '')
        except Exception as e2:
            print(f"Alternative parsing also failed: {e2}")
            print("Trying to read line by line...")
            
            # Manual parsing as last resort
            data_lines = []
            header = None
            
            with open(args.pvar, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('##'):
                        continue
                    elif line.startswith('#CHROM') or line.startswith('CHROM'):
                        header = line.replace('#', '').split('\t')
                    elif header and not line.startswith('#'):
                        data_lines.append(line.split('\t'))
            
            if not header:
                raise ValueError("Could not find proper header in PVAR file")
                
            pvar_df = pd.DataFrame(data_lines, columns=header)
    
    print(f"Found {len(pvar_df)} variants in PVAR file")
    print(f"PVAR columns: {list(pvar_df.columns)}")
    
    # Read PSAM file for population information
    print("Reading PSAM file...")
    try:
        psam_df = pd.read_csv(args.psam, sep='\t', comment='#')
        if psam_df.columns[0].startswith('#'):
            psam_df.columns = psam_df.columns.str.replace('#', '')
    except Exception as e:
        print(f"Error reading PSAM file: {e}")
        return 1
    
    print(f"PSAM columns: {list(psam_df.columns)}")
    
    # Get available populations from PSAM
    if 'SUPERPOP' in psam_df.columns:
        populations = sorted(psam_df['SUPERPOP'].unique())
        populations = [p for p in populations if p != 'UNK' and pd.notna(p)]
    else:
        print("Warning: No SUPERPOP column found, using default populations")
        populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    
    print(f"Expected populations: {populations}")
    
    # Read frequency files for each population
    freq_data = {}
    
    for freq_file in args.plink_freq:
        print(f"\nProcessing frequency file: {freq_file}")
        
        # Extract population from filename (e.g., freq_AFR.frq -> AFR)
        basename = os.path.basename(freq_file)
        if 'freq_' in basename:
            pop = basename.split('freq_')[1].split('.')[0]  # Extract pop from freq_POP.frq
        else:
            print(f"Warning: Cannot extract population from filename {freq_file}")
            continue
            
        print(f"  Population: {pop}")
        
        if not os.path.exists(freq_file):
            print(f"  Warning: File does not exist: {freq_file}")
            freq_data[pop] = {}
            continue
        
        try:
            # Read PLINK .frq file
            # Format: CHR SNP A1 A2 MAF NCHROBS
            freq_df = pd.read_csv(freq_file, sep=r'\s+')
            
            print(f"  Found {len(freq_df)} variants for {pop}")
            print(f"  Columns: {list(freq_df.columns)}")
            
            # Store frequency data indexed by SNP ID
            # Use MAF column (Minor Allele Frequency)
            if 'MAF' in freq_df.columns:
                freq_data[pop] = freq_df.set_index('SNP')['MAF'].to_dict()
            elif 'FREQ' in freq_df.columns:
                freq_data[pop] = freq_df.set_index('SNP')['FREQ'].to_dict()
            else:
                print(f"  Warning: No frequency column found in {freq_file}")
                print(f"  Available columns: {list(freq_df.columns)}")
                freq_data[pop] = {}
                
        except Exception as e:
            print(f"  Error reading {freq_file}: {e}")
            freq_data[pop] = {}
    
    # Ensure we have data for all expected populations
    for pop in populations:
        if pop not in freq_data:
            print(f"Warning: No frequency data found for population {pop}")
            freq_data[pop] = {}
    
    print(f"\nFrequency data loaded for populations: {list(freq_data.keys())}")
    
    # Create iAdmix format output
    print("Creating iAdmix format file...")
    
    # Header: #chrom position rsid A1 A2 + populations
    header_cols = ['#chrom', 'position', 'rsid', 'A1', 'A2'] + populations
    
    output_data = []
    processed_variants = 0
    skipped_variants = 0
    
    for idx, row in pvar_df.iterrows():
        if idx % 10000 == 0:
            print(f"Processed {idx} variants...")
        
        variant_id = row['ID']
        chrom = str(row['CHROM']).replace('chr', '')  # Remove chr prefix for iAdmix
        pos = row['POS']
        ref_allele = row['REF']
        alt_allele = row['ALT']
        
        # Build output row
        out_row = [chrom, pos, variant_id, ref_allele, alt_allele]
        
        # Add population frequencies
        has_freq_data = False
        for pop in populations:
            if pop in freq_data and variant_id in freq_data[pop]:
                freq = freq_data[pop][variant_id]
                # Ensure frequency is a valid number
                if pd.notna(freq) and isinstance(freq, (int, float)):
                    out_row.append(f"{float(freq):.5f}")
                    has_freq_data = True
                else:
                    out_row.append("0.00000")
            else:
                out_row.append("0.00000")  # Default frequency if not found
        
        # Only include variants that have frequency data for at least one population
        if has_freq_data:
            output_data.append(out_row)
            processed_variants += 1
        else:
            skipped_variants += 1
    
    # Write output file
    print(f"\nWriting iAdmix frequency file to {args.output}")
    
    with open(args.output, 'w') as f:
        # Write header
        f.write('\t'.join(header_cols) + '\n')
        
        # Write data
        for row in output_data:
            f.write('\t'.join(map(str, row)) + '\n')
    
    print("iAdmix frequency file created successfully!")
    print(f"Populations included: {populations}")
    print(f"Total variants processed: {processed_variants}")
    print(f"Variants skipped (no frequency data): {skipped_variants}")
    
    # Summary statistics
    print("\nFrequency data availability by population:")
    for pop in populations:
        if pop in freq_data:
            available = len([row for row in output_data if len(row) > header_cols.index(pop) and float(row[header_cols.index(pop)]) > 0])
            print(f"  {pop}: {available} variants ({available/len(output_data)*100:.1f}% if data exists)")
        else:
            print(f"  {pop}: No data")

if __name__ == '__main__':
    main()