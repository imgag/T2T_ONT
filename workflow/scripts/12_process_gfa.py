#!/usr/bin/env python3
"""
Process GFA file to update node labels according to scfmap and add color tags.
"""

import sys
import csv
import argparse
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Process GFA file with scfmap and colors')
    parser.add_argument('-i', '--gfa', required=True, help='Input GFA file')
    parser.add_argument('-s', '--scfmap', required=True, help='Scaffold map file')
    parser.add_argument('-c', '--colors', required=True, help='Colors CSV file')
    parser.add_argument('-o', '--output', required=True, help='Output GFA file')
    return parser.parse_args()

def load_scfmap(scfmap_file):
    """Load the scaffold map file and create a mapping from piece IDs to contig names."""
    id_to_contig = {}
    old_id = None
    
    with open(scfmap_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('path'):
                # Extract contig name (e.g., "contig-0000001" from "path contig-0000001 utig4-1")
                parts = line.split()
                if len(parts) >= 2:
                    old_id = parts[2]
                    contig_id = parts[1]    
                    id_to_contig[old_id] = contig_id
                
    return id_to_contig

def load_colors(colors_file):
    """Load colors from CSV file."""
    contig_colors = {}
    
    with open(colors_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('contig\tcolor'):
                continue
            parts = line.split()
            contig_id = parts[0]
            color = parts[1]
            chr = parts[2]
            hap = parts[3]
            contig_colors[contig_id] = (color, chr, hap)
                
    return contig_colors

def process_gfa(gfa_file, scfmap, colors, output_file):
    """Process GFA file to update node labels and add color tags."""
    with open(gfa_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('S'):  # Segment line
                parts = line.strip().split('\t')
                node_id = parts[1]
                
                # Check if this node ID is in the scfmap
                if node_id in scfmap:
                    # Replace node ID with contig name
                    contig_name = scfmap[node_id]
                    parts[1] = contig_name
                    
                    # Add color tag if available
                    if contig_name in colors:
                        (color, chr, hap) = colors[contig_name]
                        parts.append(f"CB:Z:{color}")
                        parts.append(f"CH:Z:{chr}")
                        parts.append(f"HP:Z:{hap}")
                
                # Write modified line
                outfile.write('\t'.join(parts) + '\n')
            else:
                # For non-segment lines, update any references to node IDs
                if line.startswith('L'):  # Link line
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        # Update source node ID
                        if parts[1] in scfmap:
                            parts[1] = scfmap[parts[1]]
                        # Update target node ID
                        if parts[3] in scfmap:
                            parts[3] = scfmap[parts[3]]
                    outfile.write('\t'.join(parts) + '\n')
                else:
                    # Pass through other lines unchanged
                    outfile.write(line)

def main():
    args = parse_args()
    
    # Load mapping from scfmap file
    scfmap = load_scfmap(args.scfmap)
    
    # Load colors from CSV file
    colors = load_colors(args.colors)
    
    # Process the GFA file
    process_gfa(args.gfa, scfmap, colors, args.output)
    
    print(f"Processed GFA file saved to {args.output}")

if __name__ == "__main__":
    main()
