#!/usr/bin/env python3
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/24_update_1000G.py

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Update 1000G PSAM file with population information')
    parser.add_argument('--psam', required=True, help='Input PSAM file')  # UPDATED: --psam instead of --fam
    parser.add_argument('--metadata', required=True, help='1000G metadata file')  
    parser.add_argument('--output', required=True, help='Output PSAM file')  # UPDATED: output PSAM
    
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
    
    # Add missing columns with default values if they don't exist
    if 'PAT' not in psam.columns:
        psam['PAT'] = '0'
    if 'MAT' not in psam.columns:
        psam['MAT'] = '0'
    if 'SEX' not in psam.columns:
        psam['SEX'] = '0'
    
    # Read 1000G metadata
    try:
        metadata = pd.read_csv(args.metadata, sep='\t')
    except Exception as e:
        print(f"Error reading metadata file: {e}")
        return
    
    # Check if required columns exist in metadata
    if 'sample' not in metadata.columns:
        print("Error: 'sample' column not found in metadata file")
        return
    
    # Merge on sample ID
    psam_updated = psam.merge(metadata[['sample', 'pop', 'super_pop']], 
                             left_on='IID', right_on='sample', how='left')
    
    # Map super population to numeric codes for phenotype
    pop_mapping = {
        'AFR': 1,  # African
        'AMR': 2,  # Admixed American  
        'EAS': 3,  # East Asian
        'EUR': 4,  # European
        'SAS': 5   # South Asian
    }
    
    # Update or create PHENO1 column with super population codes
    psam_updated['PHENO1'] = psam_updated['super_pop'].map(pop_mapping).fillna(-9).astype(int)
    
    # Also add population information as additional columns
    psam_updated['POP'] = psam_updated['pop'].fillna('UNK')
    psam_updated['SUPERPOP'] = psam_updated['super_pop'].fillna('UNK')
    
    # Select and order columns for output PSAM file
    output_cols = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO1', 'POP', 'SUPERPOP']
    
    # Ensure all output columns exist
    for col in output_cols:
        if col not in psam_updated.columns:
            if col == 'PHENO1':
                psam_updated[col] = -9
            elif col in ['POP', 'SUPERPOP']:
                psam_updated[col] = 'UNK'
            else:
                psam_updated[col] = '0'
    
    # Create final output dataframe
    output_psam = psam_updated[output_cols].copy()
    
    # Rename FID column to have # prefix (PSAM standard)
    output_psam = output_psam.rename(columns={'FID': '#FID'})
    
    # Save updated PSAM file
    output_psam.to_csv(args.output, sep='\t', index=False)
    
    # Print summary statistics
    print(f"Updated PSAM file saved to {args.output}")
    print(f"Total samples: {len(output_psam)}")
    print("Super population distribution:")
    superpop_counts = output_psam['SUPERPOP'].value_counts()
    for pop, count in superpop_counts.items():
        print(f"  {pop}: {count}")
    
    # Check for samples without population assignment
    unknown_count = (output_psam['PHENO1'] == -9).sum()
    if unknown_count > 0:
        print(f"Warning: {unknown_count} samples could not be assigned to a super population")

if __name__ == '__main__':
    main()