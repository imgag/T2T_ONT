#!/usr/bin/env python3

import argparse
import pandas as pd
import gzip
import os

def read_vcf_samples(vcf_file):
    """Extract sample IDs from VCF file"""
    print(f"Reading VCF file: {vcf_file}")
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Header line with sample names
                fields = line.strip().split('\t')
                # Samples start after the first 9 fixed columns
                samples = fields[9:]
                print(f"Found {len(samples)} samples in VCF")
                return samples
    
    raise ValueError("Could not find sample header line in VCF")

def infer_population_from_sample_id(sample_id):
    """
    Infer population from 1000G sample ID format.
    1000G samples typically have format like: HG00096, NA18525, etc.
    The first 2-3 letters often indicate the population.
    """
    # Common 1000G population prefixes
    pop_prefixes = {
        'HG': 'EUR',  # Utah residents with ancestry from northern and western Europe
        'NA': 'EUR',  # European ancestry in Utah (CEU) or others
        'GM': 'AFR',  # Gambian in Western Divisions in the Gambia
        'YRI': 'AFR', # Yoruba in Ibadan, Nigeria
        'LWK': 'AFR', # Luhya in Webuye, Kenya
        'ASW': 'AFR', # Americans of African Ancestry in SW USA
        'ACB': 'AMR', # African Caribbeans in Barbados
        'MXL': 'AMR', # Mexican Ancestry from Los Angeles USA
        'PUR': 'AMR', # Puerto Ricans from Puerto Rico
        'CLM': 'AMR', # Colombians from Medellin, Colombia
        'PEL': 'AMR', # Peruvians from Lima, Peru
        'CHB': 'EAS', # Han Chinese in Beijing, China
        'JPT': 'EAS', # Japanese in Tokyo, Japan
        'CHS': 'EAS', # Southern Han Chinese
        'CDX': 'EAS', # Chinese Dai in Xishuangbanna, China
        'KHV': 'EAS', # Kinh in Ho Chi Minh City, Vietnam
        'GIH': 'SAS', # Gujarati Indian from Houston, Texas
        'PJL': 'SAS', # Punjabi from Lahore, Pakistan
        'BEB': 'SAS', # Bengali from Bangladesh
        'STU': 'SAS', # Sri Lankan Tamil from the UK
        'ITU': 'SAS', # Indian Telugu from the UK
    }
    
    # Check prefixes
    for prefix, pop in pop_prefixes.items():
        if sample_id.startswith(prefix):
            return pop
    
    # Default to UNKNOWN if can't infer
    return 'UNKNOWN'

def main():
    parser = argparse.ArgumentParser(description='Create sample mapping file for local ancestry from VCF')
    parser.add_argument('--vcf', required=True, help='Input VCF file (reference panel)')
    parser.add_argument('--output', required=True, help='Output sample map file')
    parser.add_argument('--panel', required=True, help='1000G panel file with sample metadata')
    parser.add_argument('--pop-level', choices=['population', 'superpopulation'], default='superpopulation',
                        help='Population level to use: population (e.g., GBR, FIN) or superpopulation (e.g., EUR, AFR)')
    
    args = parser.parse_args()
    
    # Get samples from VCF
    vcf_samples = read_vcf_samples(args.vcf)
    
    # Load population mapping from panel file
    print(f"Reading population mapping from: {args.panel}")
    panel_data = pd.read_csv(args.panel, sep='\t')
    
    # Display available columns
    print(f"Panel file columns: {list(panel_data.columns)}")
    
    # Determine which columns to use based on file format
    sample_col = None
    pop_col = None
    
    # Try different possible column names
    if 'Sample name' in panel_data.columns:
        sample_col = 'Sample name'
    elif 'sample' in panel_data.columns:
        sample_col = 'sample'
    else:
        raise ValueError("Could not find sample column in panel file")
    
    if args.pop_level == 'superpopulation':
        if 'Superpopulation code' in panel_data.columns:
            pop_col = 'Superpopulation code'
        elif 'super_pop' in panel_data.columns:
            pop_col = 'super_pop'
        else:
            raise ValueError("Could not find superpopulation column in panel file")
    else:  # population
        if 'Population code' in panel_data.columns:
            pop_col = 'Population code'
        elif 'pop' in panel_data.columns:
            pop_col = 'pop'
        else:
            raise ValueError("Could not find population column in panel file")
    
    print(f"Using sample column: '{sample_col}'")
    print(f"Using population column: '{pop_col}' (level: {args.pop_level})")
    
    # Create sample to population mapping
    sample_to_pop = dict(zip(panel_data[sample_col], panel_data[pop_col]))
    print(f"Loaded {len(sample_to_pop)} sample-to-population mappings")
    
    # Create sample map dataframe matching VCF samples
    sample_data = []
    missing_samples = []
    
    for sample in vcf_samples:
        if sample in sample_to_pop:
            pop = sample_to_pop[sample]
        else:
            pop = 'UNKNOWN'
            missing_samples.append(sample)
        sample_data.append({'sample_id': sample, 'population': pop})
    
    if missing_samples:
        print(f"\nWarning: {len(missing_samples)} samples in VCF not found in panel file")
        print(f"First few missing: {missing_samples[:5]}")
    
    sample_map = pd.DataFrame(sample_data)
    
    # For superpopulation level, mark non-standard populations as QUERY
    if args.pop_level == 'superpopulation':
        # Standard 1000G superpopulations: AFR, AMR, EAS, EUR, SAS
        standard_pops = {'AFR', 'AMR', 'EAS', 'EUR', 'SAS'}
        sample_map.loc[~sample_map['population'].isin(standard_pops), 'population'] = 'QUERY'
    
    # Save sample map (without header for Gnomix/RFMix)
    sample_map.to_csv(args.output, sep='\t', index=False, header=False)
    
    # Also save a version with header for debugging
    debug_file = args.output.replace('.txt', '_with_header.txt')
    sample_map.to_csv(debug_file, sep='\t', index=False, header=True)
    
    print(f"\nSample map saved to: {args.output}")
    print(f"Debug file with header saved to: {debug_file}")
    print(f"Total samples: {len(sample_map)}")
    print(f"\n{args.pop_level.capitalize()} distribution:")
    pop_counts = sample_map['population'].value_counts()
    for pop, count in pop_counts.items():
        print(f"  {pop}: {count} samples")
    
    # Show first few entries for verification
    print(f"\nFirst 10 entries:")
    print(sample_map.head(10).to_string(index=False))

if __name__ == '__main__':
    main()