#!/usr/bin/env python3
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/25_create_ancestry_sample_map.py

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create sample mapping file for local ancestry')
    parser.add_argument('--psam', required=True, help='Input PSAM file')  # UPDATED: --psam instead of --fam
    parser.add_argument('--output', required=True, help='Output sample map file')
    
    args = parser.parse_args()
    
    # Read PSAM file (has header, unlike FAM)
    psam = pd.read_csv(args.psam, sep='\t')
    
    # Handle different possible column names for PSAM format
    if '#FID' in psam.columns:
        psam = psam.rename(columns={'#FID': 'FID'})
    
    # Ensure required columns exist
    required_cols = ['FID', 'IID']
    for col in required_cols:
        if col not in psam.columns:
            print(f"Error: Required column '{col}' not found in PSAM file")
            return
    
    # Handle phenotype column (could be PHENO1, PHENOTYPE, or similar)
    pheno_col = None
    for col in ['PHENO1', 'PHENOTYPE', 'PHENO']:
        if col in psam.columns:
            pheno_col = col
            break
    
    if pheno_col is None:
        print("Warning: No phenotype column found, using -9 (Unknown) for all samples")
        psam['PHENO'] = -9
    else:
        psam['PHENO'] = psam[pheno_col]
    
    # Create sample map (sample_id, population)
    # Assuming PHENO codes: 1=AFR, 2=AMR, 3=EAS, 4=EUR, 5=SAS, -9=Unknown
    pop_map = {1: 'AFR', 2: 'AMR', 3: 'EAS', 4: 'EUR', 5: 'SAS', -9: 'UNK'}
    
    sample_map = pd.DataFrame({
        'sample_id': psam['IID'],
        'population': psam['PHENO'].map(pop_map).fillna('UNK')
    })
    
    # Save sample map (no header for RFMix compatibility)
    sample_map.to_csv(args.output, sep='\t', index=False, header=False)
    
    print(f"Sample map saved to {args.output}")
    print(f"Total samples: {len(sample_map)}")
    print("Population counts:")
    print(sample_map['population'].value_counts().to_string())

if __name__ == '__main__':
    main()