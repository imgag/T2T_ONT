#!/usr/bin/env python3
"""
Summarize all 3D structure prediction results.

Aggregates metrics from all samples, tools, and regions.
"""

import argparse
import json
import sys
from pathlib import Path
import re


def extract_metadata_from_path(json_path):
    """Extract sample, chr, roi, res, tool, hp from file path."""
    path = Path(json_path)
    parts = path.parts
    
    # Expected path: analysis_other/3dstructure/{asm}/{chr}/{roi}/{res}/{tool}/{hp}/metrics.json
    metadata = {
        'sample': 'unknown',
        'chromosome': 'unknown',
        'roi': 'unknown',
        'resolution': 0,
        'tool': 'unknown',
        'haplotype': 'unknown'
    }
    
    try:
        idx = parts.index('3dstructure')
        if idx + 6 < len(parts):
            metadata['sample'] = parts[idx + 1]
            metadata['chromosome'] = parts[idx + 2]
            metadata['roi'] = parts[idx + 3]
            metadata['resolution'] = int(parts[idx + 4])
            metadata['tool'] = parts[idx + 5]
            metadata['haplotype'] = parts[idx + 6]
    except (ValueError, IndexError):
        pass
    
    return metadata


def main():
    parser = argparse.ArgumentParser(
        description='Summarize 3D structure prediction results')
    parser.add_argument('--metrics', nargs='+', required=True,
                       help='Input metrics JSON files')
    parser.add_argument('--output', required=True,
                       help='Output TSV summary file')
    
    args = parser.parse_args()
    
    # Collect all metrics
    all_results = []
    
    for metrics_file in args.metrics:
        try:
            with open(metrics_file, 'r') as f:
                metrics = json.load(f)
            
            metadata = extract_metadata_from_path(metrics_file)
            
            result = {
                **metadata,
                'n_coordinates': metrics.get('n_coordinates', 0),
                'dscc_spearman': metrics.get('dscc_spearman', 0.0),
                'dscc_pearson': metrics.get('dscc_pearson', 0.0),
                'radius_of_gyration': metrics.get('radius_of_gyration', 0.0),
                'reconstruction_loss': metrics.get('reconstruction_loss', float('inf')),
                'status': metrics.get('status', 'unknown')
            }
            
            all_results.append(result)
            
        except Exception as e:
            print(f"Warning: Could not process {metrics_file}: {e}", 
                  file=sys.stderr)
    
    # Write summary
    if all_results:
        columns = ['sample', 'chromosome', 'roi', 'resolution', 'tool', 
                  'haplotype', 'n_coordinates', 'dscc_spearman', 'dscc_pearson',
                  'radius_of_gyration', 'reconstruction_loss', 'status']
        
        with open(args.output, 'w') as f:
            f.write('\t'.join(columns) + '\n')
            
            for result in all_results:
                values = [str(result.get(col, '')) for col in columns]
                f.write('\t'.join(values) + '\n')
        
        print(f"Summary written to {args.output}")
        print(f"Total entries: {len(all_results)}")
        
        # Print summary statistics by tool
        tool_stats = {}
        for r in all_results:
            tool = r['tool']
            if tool not in tool_stats:
                tool_stats[tool] = []
            if r['dscc_spearman'] != 0.0:
                tool_stats[tool].append(r['dscc_spearman'])
        
        print("\nMean dSCC (Spearman) by tool:")
        for tool, values in sorted(tool_stats.items()):
            if values:
                mean_dscc = sum(values) / len(values)
                print(f"  {tool}: {mean_dscc:.4f} (n={len(values)})")
    else:
        print("Warning: No valid metrics found", file=sys.stderr)
        with open(args.output, 'w') as f:
            f.write("# No valid metrics found\n")


if __name__ == '__main__':
    main()
