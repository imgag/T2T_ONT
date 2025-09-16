#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Create sample mapping file for local ancestry')
    parser.add_argument('--fam', required=True, help='Input FAM file')
    parser.add_argument('--output', required=True, help='Output sample map file')
    
    args = parser.parse_args()
    
    # Read FAM file
    fam = pd.read_csv(args.fam, sep='\t', header=None,
                     names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    
    # Create sample map (sample_id, population)
    # Assuming PHENO codes: 1=AFR, 2=AMR, 3=EAS, 4=EUR, 5=SAS, -9=Unknown
    pop_map = {1: 'AFR', 2: 'AMR', 3: 'EAS', 4: 'EUR', 5: 'SAS', -9: 'UNK'}
    
    sample_map = pd.DataFrame({
        'sample_id': fam['IID'],
        'population': fam['PHENO'].map(pop_map)
    })
    
    # Save sample map
    sample_map.to_csv(args.output, sep='\t', index=False, header=False)
    
    print(f"Sample map saved to {args.output}")

if __name__ == '__main__':
    main()