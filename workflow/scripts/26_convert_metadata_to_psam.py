#!/usr/bin/env python3
# filepath: /mnt/storage3b/projects/no_ngsd/ahthapp1_T2T_ONT/workflow/scripts/26_convert_metadata_to_psam.py

"""
Convert 1000G metadata to PSAM file, integrating panel information for population data
PRESERVES the exact sample order from the VCF file
"""
import argparse
import subprocess
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Convert 1000G metadata to PSAM format with panel info')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--metadata', required=True, help='Input metadata file (PED format)')
    parser.add_argument('--panel', required=True, help='Input panel file with population info')
    parser.add_argument('--psam', required=True, help='Output PSAM file')
    
    args = parser.parse_args()
    
    # Get VCF sample IDs using bcftools - PRESERVE ORDER
    print("Extracting sample IDs from VCF (preserving order)...", file=sys.stderr)
    try:
        result = subprocess.run(['bcftools', 'query', '-l', args.vcf], 
                               capture_output=True, text=True, check=True)
        vcf_samples_ordered = result.stdout.strip().split('\n')  # LIST, not set!
        vcf_samples_set = set(vcf_samples_ordered)  # For fast lookup
        print(f"Found {len(vcf_samples_ordered)} samples in VCF", file=sys.stderr)
        
        # Check for duplicates in VCF
        if len(vcf_samples_ordered) != len(vcf_samples_set):
            print(f"WARNING: VCF has duplicate samples! Unique: {len(vcf_samples_set)}", file=sys.stderr)
            
    except subprocess.CalledProcessError as e:
        print(f"Error running bcftools: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Read panel file for population information
    print("Reading panel file...", file=sys.stderr)
    try:
        panel_df = pd.read_csv(args.panel, sep='\t', header=0)
        print(f"Found {len(panel_df)} samples in panel file", file=sys.stderr)
        print(f"Panel columns: {list(panel_df.columns)}", file=sys.stderr)
        
        # Rename columns to standard names if needed
        panel_df = panel_df.rename(columns={
            'sample': 'sample_id',
            'pop': 'population', 
            'super_pop': 'superpopulation',
            'gender': 'sex'
        })
        
        # Convert to dictionary for fast lookup
        panel_dict = {}
        for _, row in panel_df.iterrows():
            panel_dict[row['sample_id']] = {
                'population': row['population'],
                'superpopulation': row['superpopulation'], 
                'sex': row['sex']
            }
        
        print(f"Panel data indexed for {len(panel_dict)} samples", file=sys.stderr)
        
    except Exception as e:
        print(f"Error reading panel file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Read 1000G metadata file (PED format) - may be optional if panel has all needed info
    metadata_df = None
    metadata_dict = {}
    if args.metadata:
        print("Reading 1000G metadata file...", file=sys.stderr)
        try:
            metadata_df = pd.read_csv(args.metadata, sep='\t', header=0)
            print(f"Found {len(metadata_df)} samples in metadata", file=sys.stderr)
            print(f"Metadata columns: {list(metadata_df.columns)}", file=sys.stderr)
            
            # Create lookup dictionary for metadata (assume IID is column 1)
            for _, row in metadata_df.iterrows():
                sample_id = str(row.iloc[1])  # IID column
                metadata_dict[sample_id] = {
                    'fid': str(row.iloc[0]) if pd.notna(row.iloc[0]) else sample_id,
                    'pat_id': str(row.iloc[2]) if pd.notna(row.iloc[2]) else '0',
                    'mat_id': str(row.iloc[3]) if pd.notna(row.iloc[3]) else '0',
                    'pheno': int(row.iloc[5]) if pd.notna(row.iloc[5]) else -9
                }
                
        except Exception as e:
            print(f"Warning: Error reading metadata file: {e}", file=sys.stderr)
            print("Continuing with panel data only...", file=sys.stderr)
    
    # Create PSAM data IN THE SAME ORDER as VCF samples
    print("Creating PSAM data in VCF sample order...", file=sys.stderr)
    psam_data = []
    found_samples = []
    missing_samples = []
    
    # Process samples in VCF order
    for sample_id in vcf_samples_ordered:
        if sample_id in panel_dict:
            panel_info = panel_dict[sample_id]
            
            # Convert sex coding: male->1, female->2, unknown->0
            sex_map = {'male': 1, 'female': 2, 'M': 1, 'F': 2}
            sex_code = sex_map.get(panel_info['sex'], 0)
            
            # Get additional info from metadata if available
            if sample_id in metadata_dict:
                meta_info = metadata_dict[sample_id]
                fid = meta_info['fid']
                pat_id = meta_info['pat_id']
                mat_id = meta_info['mat_id']
                pheno = meta_info['pheno']
            else:
                # Default values if no metadata
                fid = sample_id
                pat_id = '0'
                mat_id = '0'
                pheno = -9
            
            psam_data.append({
                'FID': fid,
                'IID': sample_id,
                'PAT': pat_id,
                'MAT': mat_id,
                'SEX': sex_code,
                'PHENO1': pheno,
                'POP': panel_info['population'],
                'SUPERPOP': panel_info['superpopulation']
            })
            found_samples.append(sample_id)
        else:
            missing_samples.append(sample_id)
            print(f"Warning: Sample {sample_id} not found in panel file", file=sys.stderr)
    
    print(f"Successfully matched {len(found_samples)} samples", file=sys.stderr)
    print(f"Missing from panel: {len(missing_samples)} samples", file=sys.stderr)
    
    if len(missing_samples) > 0 and len(missing_samples) <= 10:
        print(f"Missing samples: {missing_samples}", file=sys.stderr)
    
    # Verify order preservation
    psam_sample_order = [row['IID'] for row in psam_data]
    vcf_found_order = [s for s in vcf_samples_ordered if s in panel_dict]
    
    if psam_sample_order == vcf_found_order:
        print("✅ Sample order preserved correctly", file=sys.stderr)
    else:
        print("❌ ERROR: Sample order not preserved!", file=sys.stderr)
        print(f"Expected: {vcf_found_order[:5]}...", file=sys.stderr)
        print(f"Got: {psam_sample_order[:5]}...", file=sys.stderr)
    
    # Create population summaries
    if len(psam_data) > 0:
        psam_df = pd.DataFrame(psam_data)
        
        pop_counts = psam_df['POP'].value_counts()
        print(f"Population distribution:", file=sys.stderr)
        for pop, count in pop_counts.items():
            print(f"  {pop}: {count}", file=sys.stderr)
        
        superpop_counts = psam_df['SUPERPOP'].value_counts()
        print(f"Super population distribution:", file=sys.stderr)
        for superpop, count in superpop_counts.items():
            print(f"  {superpop}: {count}", file=sys.stderr)
        
        sex_counts = psam_df['SEX'].value_counts()
        print(f"Sex distribution:", file=sys.stderr)
        for sex, count in sex_counts.items():
            sex_label = {1: 'Male', 2: 'Female', 0: 'Unknown'}.get(sex, str(sex))
            print(f"  {sex_label}: {count}", file=sys.stderr)
    
    # Write PSAM file - samples will be in VCF order
    print("Writing PSAM file in VCF sample order...", file=sys.stderr)
    header_cols = ['#FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO1', 'POP', 'SUPERPOP']
    
    with open(args.psam, 'w') as f:
        # Write header
        f.write('\t'.join(header_cols) + '\n')
        
        # Write data in order
        for row_data in psam_data:
            line_data = [
                row_data['FID'],
                row_data['IID'], 
                row_data['PAT'],
                row_data['MAT'],
                str(row_data['SEX']),
                str(row_data['PHENO1']),
                row_data['POP'],
                row_data['SUPERPOP']
            ]
            f.write('\t'.join(line_data) + '\n')
    
    # Report final statistics
    print(f"\nFinal PSAM file statistics:", file=sys.stderr)
    print(f"VCF samples: {len(vcf_samples_ordered)}", file=sys.stderr)
    print(f"Panel samples: {len(panel_dict)}", file=sys.stderr)
    if metadata_df is not None:
        print(f"Metadata samples: {len(metadata_df)}", file=sys.stderr)
    print(f"Final PSAM samples: {len(psam_data)}", file=sys.stderr)
    print(f"PSAM columns: {header_cols}", file=sys.stderr)
    
    # Report coverage
    panel_samples = set(panel_dict.keys())
    matched_samples = vcf_samples_set & panel_samples
    vcf_only = vcf_samples_set - panel_samples
    panel_only = panel_samples - vcf_samples_set
    
    print(f"\nSample overlap analysis:", file=sys.stderr)
    print(f"Samples in both VCF and panel: {len(matched_samples)}", file=sys.stderr)
    print(f"Samples in VCF but not panel: {len(vcf_only)}", file=sys.stderr)
    print(f"Samples in panel but not VCF: {len(panel_only)}", file=sys.stderr)
    
    if len(vcf_only) > 0 and len(vcf_only) <= 10:
        print(f"VCF-only samples: {list(vcf_only)}", file=sys.stderr)
    
    if len(psam_data) == 0:
        print("ERROR: No samples matched between VCF and panel!", file=sys.stderr)
        sys.exit(1)
    
    # Final verification
    print(f"\n✅ PSAM file created successfully with {len(psam_data)} samples in VCF order!", file=sys.stderr)
    print(f"First 5 samples: {[row['IID'] for row in psam_data[:5]]}", file=sys.stderr)
    print(f"Last 5 samples: {[row['IID'] for row in psam_data[-5:]]}", file=sys.stderr)

if __name__ == '__main__':
    main()