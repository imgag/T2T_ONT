#!/usr/bin/env python3

import argparse
import pandas as pd
import re

def parse_pedigree_from_names(sample_names):
    """
    Parse pedigree relationships from sample names.
    Convention: T2T04 (child), T2T04_1 (father), T2T04_2 (mother)
    """
    pedigree = {}
    
    # First, identify all base family IDs and their members
    families = {}
    
    for sample in sample_names:
        # Check if sample follows the pattern: baseID or baseID_1 or baseID_2
        if '_' in sample:
            base_id, suffix = sample.rsplit('_', 1)
            if suffix in ['1', '2']:  # Parent
                if base_id not in families:
                    families[base_id] = {'child': base_id, 'parents': {}}
                families[base_id]['parents'][suffix] = sample
        else:
            # This could be a child (if parents exist) or unrelated individual
            if sample not in families:
                families[sample] = {'child': sample, 'parents': {}}
    
    # Build pedigree structure
    for family_id, family_info in families.items():
        child = family_info['child']
        parents = family_info['parents']
        
        # Set up child
        father = parents.get('1', '0')  # '1' suffix = father
        mother = parents.get('2', '0')  # '2' suffix = mother
        
        pedigree[child] = {
            'FID': family_id,
            'IID': child,
            'PAT': father if father in sample_names else '0',
            'MAT': mother if mother in sample_names else '0'
        }
        
        # Set up parents (if they exist)
        for parent_suffix, parent_id in parents.items():
            if parent_id in sample_names:
                pedigree[parent_id] = {
                    'FID': family_id,
                    'IID': parent_id,
                    'PAT': '0',
                    'MAT': '0'
                }
    
    return pedigree

def update_fam_with_pedigree(fam_df, pedigree_info):
    """
    Update FAM file with pedigree information
    """
    updated_fam = fam_df.copy()
    
    for idx, row in updated_fam.iterrows():
        sample_id = row['IID']
        if sample_id in pedigree_info:
            # Update family ID and parental information
            updated_fam.at[idx, 'FID'] = pedigree_info[sample_id]['FID']
            updated_fam.at[idx, 'PAT'] = pedigree_info[sample_id]['PAT']
            updated_fam.at[idx, 'MAT'] = pedigree_info[sample_id]['MAT']
    
    return updated_fam

def main():
    parser = argparse.ArgumentParser(description='Update FAM file with pedigree information')
    parser.add_argument('--fam', required=True, help='Input FAM file')
    parser.add_argument('--output', required=True, help='Output FAM file')
    
    args = parser.parse_args()
    
    # Read FAM file
    fam_columns = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO']
    fam_df = pd.read_csv(args.fam, sep='\s+', header=None, names=fam_columns)
    
    print(f"Read {len(fam_df)} samples from FAM file")
    
    # Extract sample names
    sample_names = fam_df['IID'].tolist()
    
    # Parse pedigree information
    pedigree_info = parse_pedigree_from_names(sample_names)
    
    # Report pedigree structure
    families = {}
    for sample_id, info in pedigree_info.items():
        fid = info['FID']
        if fid not in families:
            families[fid] = []
        families[fid].append(sample_id)
    
    print(f"Identified {len(families)} families:")
    for fid, members in families.items():
        print(f"  Family {fid}: {', '.join(members)}")
        
        # Show parent-child relationships
        for member in members:
            if pedigree_info[member]['PAT'] != '0' or pedigree_info[member]['MAT'] != '0':
                pat = pedigree_info[member]['PAT']
                mat = pedigree_info[member]['MAT']
                print(f"    {member} -> Father: {pat}, Mother: {mat}")
    
    # Update FAM file
    updated_fam = update_fam_with_pedigree(fam_df, pedigree_info)
    
    # Save updated FAM file
    updated_fam.to_csv(args.output, sep='\t', header=False, index=False)
    
    print(f"Updated FAM file saved to {args.output}")
    
    # Summary statistics
    trios = sum(1 for info in pedigree_info.values() 
                if info['PAT'] != '0' and info['MAT'] != '0')
    parent_child_pairs = sum(1 for info in pedigree_info.values() 
                            if info['PAT'] != '0' or info['MAT'] != '0') - trios
    singletons = len(sample_names) - len(pedigree_info) + sum(1 for info in pedigree_info.values() 
                                                              if info['PAT'] == '0' and info['MAT'] == '0')
    
    print(f"\nPedigree summary:")
    print(f"  Complete trios: {trios}")
    print(f"  Parent-child pairs: {parent_child_pairs}")
    print(f"  Singletons: {singletons}")

if __name__ == '__main__':
    main()