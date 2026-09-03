#!/usr/bin/env python3
"""
Convert various 3D coordinate formats to PDB format.

Supports:
- XYZ format (simple x y z per line)
- TADbit output format
- Generic tab/space-separated coordinate files

This script converts to standard PDB format for visualization.
"""

import argparse
import sys
import os


def parse_xyz(input_file):
    """Parse simple XYZ format file and return coordinates."""
    coordinates = []
    
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 3:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    
                    # Optional feature value
                    feature = float(parts[3]) if len(parts) > 3 else 0.0
                    
                    coordinates.append({
                        'chrom': 'chr',
                        'pos': i * 1000,  # Dummy position
                        'x': x,
                        'y': y,
                        'z': z,
                        'feature': feature
                    })
                except ValueError:
                    continue
    
    return coordinates


def parse_generic_coords(input_file):
    """Parse generic coordinate file with optional chromosome and position columns."""
    coordinates = []
    
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            
            # Try different formats
            if len(parts) >= 5:
                # Format: chrom position x y z [feature]
                try:
                    chrom = parts[0]
                    position = int(parts[1])
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    feature = float(parts[5]) if len(parts) > 5 else 0.0
                    
                    coordinates.append({
                        'chrom': chrom,
                        'pos': position,
                        'x': x,
                        'y': y,
                        'z': z,
                        'feature': feature
                    })
                except (ValueError, IndexError):
                    continue
            elif len(parts) >= 3:
                # Format: x y z
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    
                    coordinates.append({
                        'chrom': 'chr',
                        'pos': i * 1000,
                        'x': x,
                        'y': y,
                        'z': z,
                        'feature': 0.0
                    })
                except ValueError:
                    continue
    
    return coordinates


def write_pdb(coordinates, output_file, scale=100.0, title="3D GENOME STRUCTURE"):
    """Write coordinates to PDB format."""
    
    with open(output_file, 'w') as pdb:
        # Header
        pdb.write(f"HEADER    {title}\n")
        if coordinates:
            first = coordinates[0]
            last = coordinates[-1]
            pdb.write(f"TITLE     {first['chrom']}:{first['pos']}-{last['pos']}\n")
        pdb.write(f"REMARK    Converted from coordinate file\n")
        pdb.write(f"REMARK    Coordinate scale factor: {scale}\n")
        
        # Get unique chromosomes for chain assignment
        chroms = list(dict.fromkeys([c['chrom'] for c in coordinates]))
        chrom_to_chain = {c: chr(65 + i % 26) for i, c in enumerate(chroms)}
        
        # Atoms
        atom_num = 1
        residue_num = 1
        prev_chrom = None
        
        for coord in coordinates:
            # New chain for new chromosome
            if coord['chrom'] != prev_chrom:
                if prev_chrom is not None:
                    pdb.write("TER\n")
                residue_num = 1
                prev_chrom = coord['chrom']
            
            chain = chrom_to_chain[coord['chrom']]
            x = coord['x'] * scale
            y = coord['y'] * scale
            z = coord['z'] * scale
            
            # B-factor can store the feature value
            bfactor = coord['feature'] * 100 if coord['feature'] else 0.0
            
            # Write ATOM record
            pdb.write(f"ATOM  {atom_num:5d}  CA  ALA {chain}{residue_num:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfactor:6.2f}           C\n")
            
            atom_num += 1
            residue_num += 1
        
        pdb.write("TER\n")
        
        # CONECT records for chain connectivity
        atom_num = 1
        prev_chrom = None
        
        for i, coord in enumerate(coordinates):
            if coord['chrom'] != prev_chrom:
                prev_chrom = coord['chrom']
                prev_atom = atom_num
            else:
                pdb.write(f"CONECT{prev_atom:5d}{atom_num:5d}\n")
                prev_atom = atom_num
            atom_num += 1
        
        pdb.write("END\n")


def main():
    parser = argparse.ArgumentParser(
        description='Convert coordinate files (XYZ, generic) to PDB format')
    parser.add_argument('input', help='Input coordinate file')
    parser.add_argument('output', help='Output PDB file')
    parser.add_argument('--scale', type=float, default=100.0,
                       help='Coordinate scale factor (default: 100)')
    parser.add_argument('--title', type=str, default='3D GENOME STRUCTURE',
                       help='Title for PDB header')
    parser.add_argument('--format', type=str, choices=['xyz', 'generic'], default='generic',
                       help='Input file format (default: generic)')
    
    args = parser.parse_args()
    
    try:
        # Auto-detect format or use specified
        if args.format == 'xyz' or args.input.endswith('.xyz'):
            coordinates = parse_xyz(args.input)
        else:
            coordinates = parse_generic_coords(args.input)
        
        if not coordinates:
            print(f"Warning: No coordinates found in {args.input}", file=sys.stderr)
            # Create empty PDB
            with open(args.output, 'w') as f:
                f.write(f"HEADER    {args.title} (EMPTY)\n")
                f.write("END\n")
            return
        
        write_pdb(coordinates, args.output, args.scale, args.title)
        print(f"Converted {len(coordinates)} coordinates to {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        # Create error placeholder PDB
        with open(args.output, 'w') as f:
            f.write(f"HEADER    {args.title} (ERROR)\n")
            f.write(f"REMARK    Error: {e}\n")
            f.write("END\n")
        sys.exit(1)


if __name__ == '__main__':
    main()