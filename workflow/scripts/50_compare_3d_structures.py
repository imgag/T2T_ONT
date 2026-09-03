#!/usr/bin/env python3
"""
Compare 3D structures from different methods.

Creates combined visualization and calculates similarity metrics
between structures from ASHIC, LorDG, HiC-GNN, ParticleChromo3D, and TADbit.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr


def parse_pdb(pdb_file):
    """Parse PDB file and extract CA atom coordinates."""
    coordinates = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name == 'CA':
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coordinates.append([x, y, z])
    
    return np.array(coordinates)


def extract_tool_and_hp(pdb_path):
    """Extract tool name and haplotype from path."""
    parts = Path(pdb_path).parts
    for i, part in enumerate(parts):
        if part in ['ashic', 'lordg', 'hicgnn', 'particlechromo3d', 'tadbit']:
            tool = part
            hp = parts[i+1] if i+1 < len(parts) else 'unknown'
            return tool, hp
    return 'unknown', 'unknown'


def align_structures(coords1, coords2):
    """
    Align two structures using Procrustes analysis.
    Returns aligned coords2 and RMSD.
    """
    if len(coords1) == 0 or len(coords2) == 0:
        return coords2, float('inf')
    
    # Truncate to same length if needed
    n = min(len(coords1), len(coords2))
    c1 = coords1[:n].copy()
    c2 = coords2[:n].copy()
    
    # Center both structures
    c1 -= c1.mean(axis=0)
    c2 -= c2.mean(axis=0)
    
    # Scale normalization
    s1 = np.sqrt((c1**2).sum())
    s2 = np.sqrt((c2**2).sum())
    if s1 > 0:
        c1 /= s1
    if s2 > 0:
        c2 /= s2
    
    # Optimal rotation using SVD
    H = c2.T @ c1
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Apply rotation
    c2_aligned = c2 @ R
    
    # Calculate RMSD
    rmsd = np.sqrt(((c1 - c2_aligned)**2).mean())
    
    return c2_aligned * s1 + coords1[:n].mean(axis=0), rmsd


def calculate_distance_correlation(coords1, coords2):
    """
    Calculate Spearman correlation of pairwise distance matrices.
    This is a measure of structural similarity independent of alignment.
    """
    if len(coords1) < 3 or len(coords2) < 3:
        return 0.0
    
    # Truncate to same length
    n = min(len(coords1), len(coords2))
    c1 = coords1[:n]
    c2 = coords2[:n]
    
    # Calculate distance matrices
    dist1 = cdist(c1, c1)
    dist2 = cdist(c2, c2)
    
    # Get upper triangular (excluding diagonal)
    idx = np.triu_indices(n, k=1)
    d1_flat = dist1[idx]
    d2_flat = dist2[idx]
    
    # Spearman correlation
    corr, _ = spearmanr(d1_flat, d2_flat)
    
    return corr if not np.isnan(corr) else 0.0


def create_combined_visualization(structures, title, output_html):
    """Create combined interactive visualization."""
    
    # Create figure with subplots
    n_structures = len(structures)
    
    # Create color scales for different tools
    tool_colors = {
        'ashic': 'Viridis',
        'lordg': 'Plasma',
        'hicgnn': 'Inferno',
        'particlechromo3d': 'Cividis',
        'tadbit': 'Purples'
    }
    
    hp_symbols = {
        'hp1': 'circle',
        'hp2': 'diamond'
    }
    
    fig = go.Figure()
    
    for i, (name, coords) in enumerate(structures.items()):
        if len(coords) == 0:
            continue
        
        tool, hp = name.split('_')
        color = np.arange(len(coords))
        
        # Add line trace
        fig.add_trace(go.Scatter3d(
            x=coords[:, 0],
            y=coords[:, 1],
            z=coords[:, 2],
            mode='lines',
            line=dict(
                color=color,
                colorscale=tool_colors.get(tool, 'Rainbow'),
                width=3
            ),
            name=f'{tool} {hp}',
            legendgroup=name,
            showlegend=True
        ))
        
        # Add markers for start/end
        fig.add_trace(go.Scatter3d(
            x=[coords[0, 0]],
            y=[coords[0, 1]],
            z=[coords[0, 2]],
            mode='markers',
            marker=dict(
                size=8, 
                color='green', 
                symbol=hp_symbols.get(hp, 'circle')
            ),
            name=f'{name} start',
            legendgroup=name,
            showlegend=False
        ))
    
    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=18)),
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='cube'
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    
    fig.write_html(output_html)


def main():
    parser = argparse.ArgumentParser(
        description='Compare 3D structures from different methods')
    parser.add_argument('--pdbs', nargs='+', required=True,
                       help='Input PDB files')
    parser.add_argument('--output-html', required=True,
                       help='Output HTML visualization')
    parser.add_argument('--output-comparison', required=True,
                       help='Output comparison TSV file')
    parser.add_argument('--title', default='3D Structure Comparison',
                       help='Title for visualization')
    
    args = parser.parse_args()
    
    # Load all structures
    structures = {}
    for pdb_file in args.pdbs:
        try:
            coords = parse_pdb(pdb_file)
            tool, hp = extract_tool_and_hp(pdb_file)
            name = f"{tool}_{hp}"
            structures[name] = coords
            print(f"Loaded {name}: {len(coords)} coordinates")
        except Exception as e:
            print(f"Warning: Could not load {pdb_file}: {e}", file=sys.stderr)
            tool, hp = extract_tool_and_hp(pdb_file)
            structures[f"{tool}_{hp}"] = np.array([])
    
    # Calculate pairwise comparisons
    comparisons = []
    structure_names = list(structures.keys())
    
    for i, name1 in enumerate(structure_names):
        for name2 in structure_names[i+1:]:
            coords1 = structures[name1]
            coords2 = structures[name2]
            
            if len(coords1) > 0 and len(coords2) > 0:
                _, rmsd = align_structures(coords1, coords2)
                dist_corr = calculate_distance_correlation(coords1, coords2)
            else:
                rmsd = float('nan')
                dist_corr = float('nan')
            
            comparisons.append({
                'structure1': name1,
                'structure2': name2,
                'rmsd': rmsd,
                'distance_correlation': dist_corr
            })
    
    # Write comparison results
    with open(args.output_comparison, 'w') as f:
        f.write("structure1\tstructure2\trmsd\tdistance_correlation\n")
        for comp in comparisons:
            f.write(f"{comp['structure1']}\t{comp['structure2']}\t"
                   f"{comp['rmsd']:.4f}\t{comp['distance_correlation']:.4f}\n")
    
    print(f"Wrote comparison to {args.output_comparison}")
    
    # Create visualization
    create_combined_visualization(structures, args.title, args.output_html)
    print(f"Created visualization: {args.output_html}")


if __name__ == '__main__':
    main()
