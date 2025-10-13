#!/usr/bin/env python3
"""
Assembly Downloader Script
Downloads genome assemblies from S3 with the following features:
- Downloads all Verkko assemblies first
- Adds random HiFiasm assemblies to reach total of 20 (or specified number)
- Groups haplotypes by sample in folders
- Performs MD5 checksum verification
- Renames files to {sample}.haplotype{N}.fasta format and unzips them
- Downloads QC files (flagger and nucflag results)
- Allows adding more assemblies later
"""

import pandas as pd
import boto3
import os
import hashlib
import random
import argparse
import json
import gzip
import shutil
from pathlib import Path
from botocore import UNSIGNED
from botocore.config import Config
import sys

class AssemblyDownloader:
    def __init__(self, tsv_file, output_dir="assemblies", config_file="download_config.json", 
                 flagger_file=None, nucflag_file=None, qc_output_dir=None):
        self.tsv_file = tsv_file
        self.output_dir = Path(output_dir)
        self.config_file = config_file
        self.flagger_file = flagger_file
        self.nucflag_file = nucflag_file
        # Use specified QC output directory or default to output_dir/qc
        self.qc_output_dir = Path(qc_output_dir) if qc_output_dir else self.output_dir / "qc"
        self.output_dir.mkdir(exist_ok=True)
        self.qc_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize S3 client for public bucket (no credentials needed)
        self.s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
        
        # Load the data - auto-detect CSV vs TSV
        if tsv_file.endswith('.csv'):
            self.df = pd.read_csv(tsv_file)
        else:
            self.df = pd.read_csv(tsv_file, sep='\t')
        
        # Filter out CHM13 reference (haplotype 0)
        self.df = self.df[self.df['haplotype'] != 0]
        
        # Load QC metadata if provided
        self.flagger_df = None
        self.nucflag_df = None
        
        if flagger_file and os.path.exists(flagger_file):
            self.flagger_df = pd.read_csv(flagger_file)
            print(f"Loaded flagger metadata: {len(self.flagger_df)} entries")
        
        if nucflag_file and os.path.exists(nucflag_file):
            self.nucflag_df = pd.read_csv(nucflag_file)
            print(f"Loaded nucflag metadata: {len(self.nucflag_df)} entries")
        
        # Load or create config for tracking downloads
        self.config = self.load_config()
    
    def load_config(self):
        """Load download configuration or create new one"""
        if os.path.exists(self.config_file):
            with open(self.config_file, 'r') as f:
                config = json.load(f)
        else:
            config = {"downloaded_samples": []}
        
        # Ensure QC tracking keys exist
        if "downloaded_flagger" not in config:
            config["downloaded_flagger"] = []
        if "downloaded_nucflag" not in config:
            config["downloaded_nucflag"] = []
            
        return config
    
    def save_config(self):
        """Save download configuration"""
        with open(self.config_file, 'w') as f:
            json.dump(self.config, f, indent=2)
    
    def parse_s3_url(self, s3_url):
        """Parse S3 URL to extract bucket and key"""
        # Handle NaN/None values
        if pd.isna(s3_url) or s3_url is None:
            return None, None
            
        # Convert to string if it's not already
        s3_url = str(s3_url)
        
        if s3_url.startswith('s3://'):
            parts = s3_url[5:].split('/', 1)
            bucket = parts[0]
            key = parts[1] if len(parts) > 1 else ''
            return bucket, key
        return None, None
    
    def download_md5_file(self, md5_s3_url):
        """Download MD5 file and extract checksum"""
        # Handle missing MD5 URLs
        if pd.isna(md5_s3_url) or md5_s3_url is None:
            return None
            
        bucket, key = self.parse_s3_url(md5_s3_url)
        if not bucket or not key:
            return None
        
        try:
            response = self.s3_client.get_object(Bucket=bucket, Key=key)
            md5_content = response['Body'].read().decode('utf-8').strip()
            # Extract just the MD5 hash (first part before any whitespace)
            return md5_content.split()[0]
        except Exception as e:
            print(f"Error downloading MD5 file {md5_s3_url}: {e}")
            return None
    
    def calculate_md5(self, file_path):
        """Calculate MD5 checksum of a file"""
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def unzip_file(self, gz_path, fasta_path):
        """Unzip a gzipped file"""
        try:
            print(f"Unzipping {gz_path.name} to {fasta_path.name}")
            with gzip.open(gz_path, 'rb') as f_in:
                with open(fasta_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Remove the gzipped file after successful extraction
            gz_path.unlink()
            print(f"✓ Successfully unzipped and removed {gz_path.name}")
            return True
            
        except Exception as e:
            print(f"Error unzipping {gz_path}: {e}")
            return False
    
    def download_qc_file(self, s3_url, local_path):
        """Download QC file from S3 (these are typically not gzipped)"""
        bucket, key = self.parse_s3_url(s3_url)
        if not bucket or not key:
            print(f"Error: Invalid S3 URL: {s3_url}")
            return False
        
        try:
            print(f"Downloading QC file {key} to {local_path}")
            self.s3_client.download_file(bucket, key, str(local_path))
            print(f"✓ Downloaded QC file {local_path.name}")
            return True
                
        except Exception as e:
            print(f"Error downloading QC file {s3_url}: {e}")
            return False
    
    def download_file(self, s3_url, local_path, md5_s3_url):
        """Download file from S3, verify checksum, and unzip"""
        bucket, key = self.parse_s3_url(s3_url)
        if not bucket or not key:
            print(f"Error: Invalid S3 URL: {s3_url}")
            return False
        
        # Get expected MD5 checksum (handle missing MD5 URLs)
        expected_md5 = None
        if not pd.isna(md5_s3_url) and md5_s3_url is not None:
            expected_md5 = self.download_md5_file(md5_s3_url)
        
        if not expected_md5:
            print(f"Warning: Could not retrieve MD5 checksum for {local_path.name}")
        
        # Create temporary gzipped file path
        temp_gz_path = local_path.with_suffix('.fa.gz')
        
        try:
            print(f"Downloading {key} to {temp_gz_path}")
            self.s3_client.download_file(bucket, key, str(temp_gz_path))
            
            # Verify MD5 checksum if available
            if expected_md5:
                actual_md5 = self.calculate_md5(temp_gz_path)
                if actual_md5 == expected_md5:
                    print(f"✓ MD5 checksum verified for {temp_gz_path.name}")
                else:
                    print(f"✗ MD5 checksum mismatch for {temp_gz_path.name}")
                    print(f"  Expected: {expected_md5}")
                    print(f"  Actual:   {actual_md5}")
                    # Remove corrupted file
                    temp_gz_path.unlink()
                    return False
            else:
                print(f"✓ Downloaded {temp_gz_path.name} (no checksum verification)")
            
            # Unzip the file
            if self.unzip_file(temp_gz_path, local_path):
                return True
            else:
                return False
                
        except Exception as e:
            print(f"Error downloading {s3_url}: {e}")
            # Clean up any partial files
            if temp_gz_path.exists():
                temp_gz_path.unlink()
            return False
    
    def get_haplotype_name(self, haplotype, assembly_name):
        """Determine haplotype name - always use haplotype1/haplotype2 format"""
        return f'haplotype{haplotype}'
    
    def download_qc_files_for_sample(self, sample, sample_dir=None):
        """Download QC files for a specific sample to QC output directory"""
        # Create sample-specific QC directory in the specified QC output location
        sample_qc_dir = self.qc_output_dir / sample
        sample_qc_dir.mkdir(parents=True, exist_ok=True)
        
        qc_success = True
        
        # Download flagger files
        if self.flagger_df is not None:
            flagger_files = self.flagger_df[self.flagger_df['sample_id'] == sample]
            for _, row in flagger_files.iterrows():
                haplotype_name = self.get_haplotype_name(row['haplotype'], row['assembly_name'])
                filename = f"{sample}.{haplotype_name}.flagger.bed"
                local_path = sample_qc_dir / filename
                
                file_key = f"{sample}_{row['haplotype']}_flagger"
                if file_key not in self.config['downloaded_flagger']:
                    if self.download_qc_file(row['location'], local_path):
                        self.config['downloaded_flagger'].append(file_key)
                    else:
                        qc_success = False
                        print(f"Failed to download flagger file for {sample} haplotype {row['haplotype']}")
                else:
                    print(f"Skipping flagger file for {sample} haplotype {row['haplotype']} (already downloaded)")
        
        # Download nucflag files
        if self.nucflag_df is not None:
            nucflag_files = self.nucflag_df[self.nucflag_df['sample_id'] == sample]
            for _, row in nucflag_files.iterrows():
                haplotype_name = self.get_haplotype_name(row['haplotype'], row['assembly_name'])
                filename = f"{sample}.{haplotype_name}.nucflag.bed"
                local_path = sample_qc_dir / filename
                
                file_key = f"{sample}_{row['haplotype']}_nucflag"
                if file_key not in self.config['downloaded_nucflag']:
                    if self.download_qc_file(row['location'], local_path):
                        self.config['downloaded_nucflag'].append(file_key)
                    else:
                        qc_success = False
                        print(f"Failed to download nucflag file for {sample} haplotype {row['haplotype']}")
                else:
                    print(f"Skipping nucflag file for {sample} haplotype {row['haplotype']} (already downloaded)")
        
        return qc_success
    
    def select_assemblies(self, total_count=20):
        """Select assemblies to download based on criteria"""
        # Get all Verkko assemblies
        verkko_assemblies = self.df[self.df['assembly_method'] == 'verkko']
        
        # Get available HiFiasm assemblies (excluding already selected Verkko samples)
        verkko_samples = set(verkko_assemblies['sample_id'].unique())
        hifiasm_assemblies = self.df[
            (self.df['assembly_method'] == 'hifiasm') & 
            (~self.df['sample_id'].isin(verkko_samples))
        ]
        
        selected_assemblies = verkko_assemblies.copy()
        
        # Calculate how many more samples we need
        verkko_sample_count = len(verkko_samples)
        remaining_slots = total_count - verkko_sample_count
        
        if remaining_slots > 0:
            # Get available HiFiasm samples
            available_hifiasm_samples = list(hifiasm_assemblies['sample_id'].unique())
            
            # Randomly select HiFiasm samples
            if len(available_hifiasm_samples) >= remaining_slots:
                selected_hifiasm_samples = random.sample(available_hifiasm_samples, remaining_slots)
            else:
                selected_hifiasm_samples = available_hifiasm_samples
                print(f"Warning: Only {len(available_hifiasm_samples)} HiFiasm samples available, "
                      f"but {remaining_slots} were needed")
            
            # Add selected HiFiasm assemblies
            selected_hifiasm = hifiasm_assemblies[
                hifiasm_assemblies['sample_id'].isin(selected_hifiasm_samples)
            ]
            selected_assemblies = pd.concat([selected_assemblies, selected_hifiasm])
        
        return selected_assemblies
    
    def download_assemblies(self, total_count=20, force_redownload=False, qc_only=False):
        """Download selected assemblies and their QC files"""
        if qc_only:
            # For QC-only mode, only process samples that have already been downloaded
            downloaded_samples = self.config['downloaded_samples']
            if not downloaded_samples:
                print("No samples have been downloaded yet. Cannot download QC files without assemblies.")
                return
            
            # Filter the dataframe to only include already downloaded samples
            selected = self.df[self.df['sample_id'].isin(downloaded_samples)]
            
            if selected.empty:
                print("No assembly data found for downloaded samples in the current dataset.")
                return
                
            print(f"Downloading QC files for {len(selected['sample_id'].unique())} downloaded samples")
            print(f"QC files will be saved to: {self.qc_output_dir}")
        else:
            # Select assemblies to download
            selected = self.select_assemblies(total_count)
            
            print(f"Selected {len(selected['sample_id'].unique())} samples for download:")
            for assembler in selected['assembly_method'].unique():
                count = len(selected[selected['assembly_method'] == assembler]['sample_id'].unique())
                print(f"  - {assembler}: {count} samples")
        
        successful_downloads = 0
        failed_downloads = 0
        
        # Group by sample to keep haplotypes together
        for sample, group in selected.groupby('sample_id'):
            sample_dir = self.output_dir / sample
            sample_dir.mkdir(exist_ok=True)
            
            # Skip assembly downloads if already downloaded (unless force redownload) or if QC only
            if sample in self.config['downloaded_samples'] and not force_redownload and not qc_only:
                print(f"Skipping assembly download for {sample} (already downloaded)")
            elif not qc_only:
                print(f"\n--- Downloading sample: {sample} ---")
                sample_success = True
                
                for _, row in group.iterrows():
                    # Determine haplotype name
                    haplotype_name = self.get_haplotype_name(row['haplotype'], row['assembly_name'])
                    
                    # Create new filename: {sample}.haplotype{N}.fasta
                    new_filename = f"{sample}.{haplotype_name}.fasta"
                    local_path = sample_dir / new_filename
                    
                    # Download, verify, and unzip
                    if self.download_file(row['assembly'], local_path, row['assembly_md5']):
                        successful_downloads += 1
                    else:
                        failed_downloads += 1
                        sample_success = False
                
                # Mark sample as downloaded if all files succeeded
                if sample_success and sample not in self.config['downloaded_samples']:
                    self.config['downloaded_samples'].append(sample)
            
            # Download QC files ONLY if the sample is in downloaded_samples
            # This ensures we never download QC files for samples without assemblies
            if (self.flagger_df is not None or self.nucflag_df is not None) and sample in self.config['downloaded_samples']:
                print(f"\n--- Downloading QC files for sample: {sample} ---")
                self.download_qc_files_for_sample(sample, sample_dir)
            elif (self.flagger_df is not None or self.nucflag_df is not None) and sample not in self.config['downloaded_samples']:
                print(f"Skipping QC files for {sample} (assembly not downloaded)")
        
        # Save config
        self.save_config()
        
        if not qc_only:
            print(f"\n--- Download Summary ---")
            print(f"Successful downloads: {successful_downloads}")
            print(f"Failed downloads: {failed_downloads}")
            print(f"Total samples downloaded: {len(self.config['downloaded_samples'])}")
        
        if self.flagger_df is not None or self.nucflag_df is not None:
            print(f"Total flagger files downloaded: {len(self.config['downloaded_flagger'])}")
            print(f"Total nucflag files downloaded: {len(self.config['downloaded_nucflag'])}")
    
    def add_more_assemblies(self, additional_count):
        """Add more assemblies to existing collection"""
        current_total = len(self.config['downloaded_samples'])
        new_total = current_total + additional_count
        
        print(f"Adding {additional_count} more assemblies...")
        print(f"Current total: {current_total}")
        print(f"New total target: {new_total}")
        
        self.download_assemblies(new_total)
    
    def list_downloaded(self):
        """List currently downloaded samples"""
        print("Downloaded samples:")
        for i, sample in enumerate(self.config['downloaded_samples'], 1):
            print(f"{i:2d}. {sample}")
        print(f"\nTotal: {len(self.config['downloaded_samples'])} samples")
        
        if self.config['downloaded_flagger']:
            print(f"Flagger files downloaded: {len(self.config['downloaded_flagger'])}")
        if self.config['downloaded_nucflag']:
            print(f"Nucflag files downloaded: {len(self.config['downloaded_nucflag'])}")

def main():
    parser = argparse.ArgumentParser(description='Download genome assemblies and QC files from S3')
    parser.add_argument('tsv_file', help='TSV/CSV file containing assembly information')
    parser.add_argument('--output-dir', '-o', default='assemblies', 
                       help='Output directory for assemblies (default: assemblies)')
    parser.add_argument('--qc-output-dir', help='Output directory for QC files (default: output-dir/qc)')
    parser.add_argument('--total-count', '-n', type=int, default=20,
                       help='Total number of samples to download (default: 20)')
    parser.add_argument('--add-more', '-a', type=int,
                       help='Add more assemblies to existing collection')
    parser.add_argument('--force-redownload', '-f', action='store_true',
                       help='Force redownload of existing files')
    parser.add_argument('--list-downloaded', '-l', action='store_true',
                       help='List currently downloaded samples')
    parser.add_argument('--config-file', '-c', default='download_config.json',
                       help='Configuration file to track downloads')
    parser.add_argument('--flagger-file', help='CSV file with flagger QC metadata')
    parser.add_argument('--nucflag-file', help='CSV file with nucflag QC metadata')
    parser.add_argument('--qc-only', action='store_true',
                       help='Download only QC files (skip assembly downloads)')
    
    args = parser.parse_args()
    
    # Check if TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV/CSV file '{args.tsv_file}' not found")
        sys.exit(1)
    
    # Initialize downloader
    downloader = AssemblyDownloader(
        args.tsv_file, 
        args.output_dir, 
        args.config_file,
        args.flagger_file,
        args.nucflag_file,
        args.qc_output_dir
    )
    
    # Execute requested action
    if args.list_downloaded:
        downloader.list_downloaded()
    elif args.add_more:
        downloader.add_more_assemblies(args.add_more)
    else:
        downloader.download_assemblies(args.total_count, args.force_redownload, args.qc_only)

if __name__ == "__main__":
    main()