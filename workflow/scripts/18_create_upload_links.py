#!/usr/bin/env python3
import os
import yaml
from pathlib import Path

# Paths
DATASETS_YML = 'data/datasets.yml'
FINISHED_SAMPLES_YML = 'data/finished_samples.yml'
UPLOAD_ROOT = 'data/upload_folder'
BASECALLED_SUP = 'data/basecalled/sup'
BASECALLED_APK = 'data/basecalled/apk'

# Data types to process (not HQ types)
DATA_TYPES = ['UL', 'APK', 'POREC']

# Helper to find SUP or APK BAMs
def find_sup_bam(flowcell_id):
    dir_path = Path(BASECALLED_SUP) / flowcell_id
    if not dir_path.exists():
        return None
    for f in dir_path.glob('*.sup.unmapped.bam'):
        return f
    return None

def find_apk_bam(flowcell_id):
    dir_path = Path(BASECALLED_APK) / flowcell_id
    if not dir_path.exists():
        return None
    for f in dir_path.glob('*.apk.unmapped.bam'):
        return f
    return None

def main():
    # Load finished samples
    with open(FINISHED_SAMPLES_YML) as f:
        finished_samples = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    # Load datasets
    with open(DATASETS_YML) as f:
        datasets = yaml.safe_load(f)

    os.makedirs(UPLOAD_ROOT, exist_ok=True)

    for sample in finished_samples:
        if sample not in datasets:
            print(f"Warning: {sample} not found in datasets.yml")
            continue
        sample_data = datasets[sample]
        sample_upload_dir = Path(UPLOAD_ROOT) / sample
        sample_upload_dir.mkdir(parents=True, exist_ok=True)

        for dtype in DATA_TYPES:
            if dtype not in sample_data:
                continue
            dtype_upload_dir = sample_upload_dir / dtype
            dtype_upload_dir.mkdir(exist_ok=True)
            entries = sample_data[dtype]
            if not isinstance(entries, list):
                entries = [entries]
            for entry in entries:
                # Only process if entry is a raw folder (for UL) or matches a known basecalled BAM (for APK/POREC)
                # For UL and POREC, look for SUP BAM; for APK, look for APK BAM
                # Extract flowcell/sample id from path
                # e.g. data/raw/set1/24070/24070LRa002_04503 -> 24070LRa002_04503
                # e.g. data/raw/25006LRa014_05173 -> 25006LRa014_05173
                # e.g. data/raw/24070LRa_PoreC_04921 -> 24070LRa_PoreC_04921
                # e.g. data/raw/25006LRa010_04919/25006LRa010_04919 -> 25006LRa010_04919
                if dtype == 'APK':
                    # APK can be a folder or a file
                    if isinstance(entry, str):
                        # Try to extract the last part as flowcell id
                        flowcell_id = Path(entry).stem
                        bam = find_apk_bam(flowcell_id)
                        if bam:
                            link_name = dtype_upload_dir / bam.name
                            if not link_name.exists():
                                os.symlink(os.path.relpath(bam, dtype_upload_dir), link_name)
                        else:
                            print(f"Warning: APK BAM not found for {flowcell_id}")
                else:  # UL or POREC
                    if isinstance(entry, str):
                        flowcell_id = Path(entry).name
                        bam = find_sup_bam(flowcell_id)
                        if bam:
                            link_name = dtype_upload_dir / bam.name
                            if not link_name.exists():
                                os.symlink(os.path.relpath(bam, dtype_upload_dir), link_name)
                        else:
                            print(f"Warning: SUP BAM not found for {flowcell_id}")

if __name__ == '__main__':
    main() 