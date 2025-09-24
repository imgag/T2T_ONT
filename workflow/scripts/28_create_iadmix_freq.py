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
    
    # Read PVAR file for variant information
    print("Reading PVAR file...")
    try:
        # Read PVAR file by skipping comment lines (##) and finding the header
        with open(args.pvar, 'r') as f:
            lines = f.readlines()
        
        # Find the first non-comment line (header or data)
        header_line_idx = 0
        for i, line in enumerate(lines):
            if not line.strip().startswith('##'):
                header_line_idx = i
                break
        
        print(f"Found data starting at line {header_line_idx + 1}")
        first_data_line = lines[header_line_idx].strip()
        print(f"First data line: {first_data_line}")
        
        # Check if the first non-comment line is a header
        if first_data_line.startswith('#CHROM') or 'CHROM' in first_data_line:
            print("Found header line")
            # Read with header, skipping comment lines
            pvar_df = pd.read_csv(args.pvar, sep='\t', dtype=str, skiprows=header_line_idx, nrows=None, comment=None)
            # Clean column names
            pvar_df.columns = pvar_df.columns.str.replace('#', '')
        else:
            print("No header found, using standard PVAR column names")
            # Read without header, skipping comment lines
            pvar_df = pd.read_csv(args.pvar, sep='\t', header=None, dtype=str, skiprows=header_line_idx, comment=None)
            # Standard PVAR columns
            expected_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            # Assign column names based on number of columns
            if len(pvar_df.columns) >= 5:
                pvar_df.columns = expected_cols[:len(pvar_df.columns)]
            else:
                raise ValueError(f"PVAR file has only {len(pvar_df.columns)} columns, expected at least 5")
        
    except Exception as e:
        print(f"Error reading PVAR file: {e}")
        print("Trying alternative parsing method...")
        
        # Alternative: manually parse the file
        try:
            data_lines = []
            header = None
            
            with open(args.pvar, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('##') or not line:
                        continue  # Skip comment lines and empty lines
                    elif line.startswith('#CHROM') or (line.startswith('CHROM') and not line.startswith('##')):
                        header = line.replace('#', '').split('\t')
                        print(f"Found header: {header}")
                    elif not line.startswith('#'):
                        data_lines.append(line.split('\t'))
            
            if header:
                pvar_df = pd.DataFrame(data_lines, columns=header)
            else:
                # Use standard column names if no header found
                expected_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
                if data_lines and len(data_lines[0]) >= 5:
                    pvar_df = pd.DataFrame(data_lines, columns=expected_cols[:len(data_lines[0])])
                else:
                    raise ValueError("Could not parse PVAR file properly")
                    
        except Exception as e2:
            print(f"Alternative parsing also failed: {e2}")
            return 1
    
    print(f"Found {len(pvar_df)} variants in PVAR file")
    print(f"PVAR columns: {list(pvar_df.columns)}")
    
    # Verify required columns exist
    required_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
    missing_cols = [col for col in required_cols if col not in pvar_df.columns]
    if missing_cols:
        print(f"Error: Missing required columns in PVAR: {missing_cols}")
        print(f"Available columns: {list(pvar_df.columns)}")
        return 1
    
    # Read PSAM file for population information
    print("Reading PSAM file...")
    try:
        # Check the first line to see if it has a header
        with open(args.psam, 'r') as f:
            first_line = f.readline().strip()
            print(f"First line of PSAM: {first_line}")
        
        if first_line.startswith('#'):
            print("PSAM has header starting with #")
            psam_df = pd.read_csv(args.psam, sep='\t', dtype=str)
            # Clean column names by removing # prefix
            psam_df.columns = [col.lstrip('#') for col in psam_df.columns]
        else:
            print("PSAM has no # prefix in header")
            psam_df = pd.read_csv(args.psam, sep='\t', dtype=str)
            
        print(f"PSAM columns after cleaning: {list(psam_df.columns)}")
        
    except Exception as e:
        print(f"Error reading PSAM file: {e}")
        return 1
    
    # Get available populations from PSAM
    if 'SUPERPOP' in psam_df.columns:
        populations = sorted(psam_df['SUPERPOP'].unique())
        populations = [p for p in populations if p != 'UNK' and p != 'QUERY' and pd.notna(p)]
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
            pop = basename.split('freq_')[1].split('.')[0]
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
            freq_df = pd.read_csv(freq_file, sep=r'\s+', dtype={'CHR': str, 'SNP': str})
            
            print(f"  Found {len(freq_df)} variants for {pop}")
            print(f"  Columns: {list(freq_df.columns)}")
            
            # Store frequency data indexed by SNP ID
            if 'MAF' in freq_df.columns:
                # Convert MAF to float, handle any non-numeric values
                freq_df['MAF'] = pd.to_numeric(freq_df['MAF'], errors='coerce')
                freq_data[pop] = freq_df.set_index('SNP')['MAF'].to_dict()
            elif 'FREQ' in freq_df.columns:
                freq_df['FREQ'] = pd.to_numeric(freq_df['FREQ'], errors='coerce')
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
        if idx % 100000 == 0:
            print(f"Processed {idx} variants...")
        
        try:
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
                
        except Exception as e:
            print(f"Error processing variant at index {idx}: {e}")
            skipped_variants += 1
            continue
    
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