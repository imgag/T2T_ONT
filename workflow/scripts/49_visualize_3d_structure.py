#!/usr/bin/env python3
"""
Visualize 3D chromatin structure from PDB file.

Creates both static PNG and interactive HTML visualizations.
Uses matplotlib for static images and plotly for interactive.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go


def parse_pdb(pdb_file):
    """Parse PDB file and extract CA atom coordinates."""
    coordinates = []
    chains = []
    residue_nums = []
    bfactors = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Standard PDB format
                atom_name = line[12:16].strip()
                if atom_name == 'CA':  # Only CA atoms (backbone)
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    chain = line[21]
                    res_num = int(line[22:26])
                    try:
                        bfactor = float(line[60:66])
                    except:
                        bfactor = 0.0
                    
                    coordinates.append([x, y, z])
                    chains.append(chain)
                    residue_nums.append(res_num)
                    bfactors.append(bfactor)
    
    return np.array(coordinates), chains, residue_nums, bfactors


def create_colormap(n_points):
    """Create rainbow colormap for chromosome visualization."""
    colors = plt.cm.rainbow(np.linspace(0, 1, n_points))
    return colors


def create_static_visualization(coords, title, output_file):
    """Create static 3D visualization using matplotlib."""
    if len(coords) == 0:
        # Create placeholder image
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.text(0.5, 0.5, 'No coordinates found', 
                ha='center', va='center', fontsize=16)
        ax.axis('off')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        return
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get colors based on position along the chromosome
    colors = create_colormap(len(coords))
    
    # Plot the 3D structure as a line with color gradient
    for i in range(len(coords) - 1):
        ax.plot(coords[i:i+2, 0], coords[i:i+2, 1], coords[i:i+2, 2],
               color=colors[i], linewidth=2)
    
    # Add scatter points
    scatter = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                        c=np.arange(len(coords)), cmap='rainbow',
                        s=20, alpha=0.7)
    
    # Mark start and end
    ax.scatter(*coords[0], color='green', s=100, marker='^', label='Start')
    ax.scatter(*coords[-1], color='red', s=100, marker='v', label='End')
    
    # Labels and title
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title(title)
    ax.legend()
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, aspect=20)
    cbar.set_label('Genomic position')
    
    # Make equal aspect ratio
    max_range = np.array([coords[:, 0].max() - coords[:, 0].min(),
                         coords[:, 1].max() - coords[:, 1].min(),
                         coords[:, 2].max() - coords[:, 2].min()]).max() / 2.0
    mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
    mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
    mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def create_interactive_visualization(coords, chains, bfactors, title, output_file):
    """Create interactive 3D visualization using plotly."""
    if len(coords) == 0:
        # Create placeholder HTML
        fig = go.Figure()
        fig.add_annotation(text="No coordinates found",
                          xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False,
                          font=dict(size=20))
        fig.write_html(output_file)
        return
    
    # Create color scale based on position
    colors = np.arange(len(coords))
    
    # Main trace - line connecting all points
    line_trace = go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode='lines',
        line=dict(
            color=colors,
            colorscale='Rainbow',
            width=4
        ),
        name='Chromatin fiber'
    )
    
    # Scatter points with hover info
    scatter_trace = go.Scatter3d(
        x=coords[:, 0],
        y=coords[:, 1],
        z=coords[:, 2],
        mode='markers',
        marker=dict(
            size=4,
            color=colors,
            colorscale='Rainbow',
            colorbar=dict(
                title='Position',
                thickness=15
            ),
            opacity=0.8
        ),
        text=[f"Pos: {i}<br>B-factor: {b:.2f}" for i, b in enumerate(bfactors)],
        hoverinfo='text',
        name='Loci'
    )
    
    # Start marker
    start_trace = go.Scatter3d(
        x=[coords[0, 0]],
        y=[coords[0, 1]],
        z=[coords[0, 2]],
        mode='markers+text',
        marker=dict(size=10, color='green', symbol='diamond'),
        text=['Start'],
        textposition='top center',
        name='Start'
    )
    
    # End marker
    end_trace = go.Scatter3d(
        x=[coords[-1, 0]],
        y=[coords[-1, 1]],
        z=[coords[-1, 2]],
        mode='markers+text',
        marker=dict(size=10, color='red', symbol='diamond'),
        text=['End'],
        textposition='top center',
        name='End'
    )
    
    # Create figure
    fig = go.Figure(data=[line_trace, scatter_trace, start_trace, end_trace])
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=title,
            x=0.5,
            font=dict(size=18)
        ),
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            aspectmode='cube'
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        margin=dict(l=0, r=0, t=50, b=0)
    )
    
    # Write HTML
    fig.write_html(output_file)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize 3D chromatin structure from PDB file')
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--output-png', required=True, help='Output PNG file')
    parser.add_argument('--output-html', required=True, help='Output HTML file')
    parser.add_argument('--title', default='3D Chromatin Structure',
                       help='Title for the visualization')
    
    args = parser.parse_args()
    
    try:
        coords, chains, residue_nums, bfactors = parse_pdb(args.pdb)
        
        print(f"Loaded {len(coords)} coordinates from {args.pdb}")
        
        # Create visualizations
        create_static_visualization(coords, args.title, args.output_png)
        print(f"Created static visualization: {args.output_png}")
        
        create_interactive_visualization(coords, chains, bfactors, args.title, args.output_html)
        print(f"Created interactive visualization: {args.output_html}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        
        # Create placeholder outputs
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center')
        ax.axis('off')
        plt.savefig(args.output_png, dpi=150)
        plt.close()
        
        with open(args.output_html, 'w') as f:
            f.write(f"<html><body><h1>Error: {e}</h1></body></html>")
        
        sys.exit(1)


if __name__ == '__main__':
    main()
