#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Update 1000G FAM file with population information')
    parser.add_argument('--fam', required=True, help='Input FAM file')
    parser.add_argument('--metadata', required=True, help='1000G metadata file')  
    parser.add_argument('--output', required=True, help='Output FAM file')
    
    args = parser.parse_args()
    
    # Read FAM file
    fam = pd.read_csv(args.fam, sep='\t', header=None, 
                     names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    
    # Read metadata
    metadata = pd.read_csv(args.metadata, sep='\t')
    
    # Merge on sample ID
    fam_updated = fam.merge(metadata[['sample', 'pop', 'super_pop']], 
                           left_on='IID', right_on='sample', how='left')
    
    # Update phenotype column with super population
    fam_updated['PHENO'] = fam_updated['super_pop'].map({
        'AFR': 1, 'AMR': 2, 'EAS': 3, 'EUR': 4, 'SAS': 5
    }).fillna(-9)
    
    # Save updated FAM file
    fam_updated[['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO']].to_csv(
        args.output, sep='\t', header=False, index=False)
    
    print(f"Updated FAM file saved to {args.output}")

if __name__ == '__main__':
    main()