#!/usr/bin/env python3

"""
Create PSAM file from VCF and metadata
- Reference samples get POP and SUPERPOP from metadata
- Query samples get POP=QUERY and SUPERPOP=QRY
- Preserves VCF sample order
"""
import argparse
import subprocess
import pandas as pd
import sys
import gzip


def read_vcf_samples(vcf_file):
    """Extract sample IDs from VCF file preserving order"""
    print(f"Reading VCF: {vcf_file}", file=sys.stderr)
    
    # Check if gzipped
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    
    with open_func(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:]  # Samples start after column 9
                print(f"Found {len(samples)} samples in VCF", file=sys.stderr)
                return samples
    
    raise ValueError("Could not find #CHROM header line in VCF")

def main():
    parser = argparse.ArgumentParser(description='Create PSAM file from VCF and metadata')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--panel', required=True, help='Metadata file with sample population info')
    parser.add_argument('--output', required=True, help='Output PSAM file')
    
    args = parser.parse_args()
    
    # Get VCF samples in order
    vcf_samples = read_vcf_samples(args.vcf)
    
    # Read metadata file
    print(f"Reading metadata: {args.panel}", file=sys.stderr)
    metadata = pd.read_csv(args.panel, sep='\t')
    print(f"Metadata columns: {list(metadata.columns)}", file=sys.stderr)
    print(f"Metadata has {len(metadata)} samples", file=sys.stderr)
    
    # Identify column names
    sample_col = 'Sample name' if 'Sample name' in metadata.columns else 'sample'
    pop_col = 'Population code' if 'Population code' in metadata.columns else 'pop'
    superpop_col = 'Superpopulation code' if 'Superpopulation code' in metadata.columns else 'super_pop'
    sex_col = 'Sex' if 'Sex' in metadata.columns else 'gender'
    
    # Create lookup dictionary
    metadata_dict = {}
    for _, row in metadata.iterrows():
        sample_id = row[sample_col]
        sex_val = row.get(sex_col, 'unknown') if sex_col in row.index else 'unknown'
        metadata_dict[sample_id] = {
            'pop': row[pop_col],
            'superpop': row[superpop_col],
            'sex': sex_val
        }
    
    print(f"Loaded metadata for {len(metadata_dict)} samples", file=sys.stderr)
    
    # Create PSAM data in VCF order
    psam_data = []
    ref_count = 0
    query_count = 0
    
    for sample_id in vcf_samples:
        # Convert sex to PLINK format: 1=male, 2=female, 0=unknown
        if sample_id in metadata_dict:
            # Reference sample
            info = metadata_dict[sample_id]
            sex = info['sex'].lower()
            sex_code = 1 if sex in ['male', 'm'] else (2 if sex in ['female', 'f'] else 0)
            
            psam_data.append({
                'FID': sample_id,
                'IID': sample_id,
                'PAT': '0',
                'MAT': '0',
                'SEX': sex_code,
                'POP': info['pop'],
                'SUPERPOP': info['superpop']
            })
            ref_count += 1
        else:
            # Query sample (our T2T samples)
            psam_data.append({
                'FID': sample_id,
                'IID': sample_id,
                'PAT': '0',
                'MAT': '0',
                'SEX': 0,  # Unknown sex for query samples
                'POP': 'QUERY',
                'SUPERPOP': 'QRY'
            })
            query_count += 1
    
    print(f"\nSample classification:", file=sys.stderr)
    print(f"  Reference samples: {ref_count}", file=sys.stderr)
    print(f"  Query samples: {query_count}", file=sys.stderr)
    print(f"  Total samples: {len(psam_data)}", file=sys.stderr)
    
    # Create DataFrame
    psam_df = pd.DataFrame(psam_data)
    
    # Show distribution
    print(f"\nSuperpopulation distribution:", file=sys.stderr)
    for superpop, count in psam_df['SUPERPOP'].value_counts().items():
        print(f"  {superpop}: {count}", file=sys.stderr)
    
    # Write PSAM file with proper header
    print(f"\nWriting PSAM file: {args.output}", file=sys.stderr)
    with open(args.output, 'w') as f:
        # Write header with # prefix
        f.write('#FID\tIID\tPAT\tMAT\tSEX\tPOP\tSUPERPOP\n')
        
        # Write data
        for row in psam_data:
            f.write(f"{row['FID']}\t{row['IID']}\t{row['PAT']}\t{row['MAT']}\t{row['SEX']}\t{row['POP']}\t{row['SUPERPOP']}\n")
    
    print(f"First 5 samples: {[r['IID'] for r in psam_data[:5]]}", file=sys.stderr)
    print(f"Last 5 samples: {[r['IID'] for r in psam_data[-5:]]}", file=sys.stderr)

if __name__ == '__main__':
    main()