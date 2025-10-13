#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create sample mapping file for local ancestry')
    parser.add_argument('--psam', required=True, help='Input PSAM file')
    parser.add_argument('--output', required=True, help='Output sample map file')
    
    args = parser.parse_args()
    
    print(f"Reading PSAM file: {args.psam}")
    
    # Read PSAM file
    psam = pd.read_csv(args.psam, sep='\t')
    
    # Clean column names (remove # prefix if present)
    psam.columns = [col.lstrip('#') for col in psam.columns]
    
    print(f"PSAM columns: {list(psam.columns)}")
    print(f"Total samples: {len(psam)}")
    
    # Check if required columns exist
    if 'IID' not in psam.columns:
        print("Error: IID column not found in PSAM file")
        return 1
    
    # Determine population column
    if 'SUPERPOP' in psam.columns:
        pop_col = 'SUPERPOP'
    elif 'POP' in psam.columns:
        pop_col = 'POP'
    else:
        print("Error: No population column (SUPERPOP or POP) found in PSAM file")
        return 1
    
    print(f"Using population column: {pop_col}")
    
    # Create sample map
    sample_map = pd.DataFrame({
        'sample_id': psam['IID'],
        'population': psam[pop_col].fillna('UNK')
    })
    
    # Map query samples (samples that don't have standard 1000G population codes)
    # Standard 1000G superpopulations: AFR, AMR, EAS, EUR, SAS
    standard_pops = {'AFR', 'AMR', 'EAS', 'EUR', 'SAS'}
    
    # Mark samples not in standard populations as QUERY
    sample_map.loc[~sample_map['population'].isin(standard_pops), 'population'] = 'QUERY'
    
    # Save sample map (with header for debugging, but RFMix might need without header)
    sample_map.to_csv(args.output, sep='\t', index=False, header=False)
    
    # Also save a version with header for debugging
    debug_file = args.output.replace('.txt', '_with_header.txt')
    sample_map.to_csv(debug_file, sep='\t', index=False, header=True)
    
    print(f"Sample map saved to: {args.output}")
    print(f"Debug file with header saved to: {debug_file}")
    print(f"Total samples: {len(sample_map)}")
    print("\nPopulation distribution:")
    pop_counts = sample_map['population'].value_counts()
    for pop, count in pop_counts.items():
        print(f"  {pop}: {count} samples")
    
    # Show first few entries for verification
    print(f"\nFirst 10 entries:")
    print(sample_map.head(10).to_string(index=False))

if __name__ == '__main__':
    main()