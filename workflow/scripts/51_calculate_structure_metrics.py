#!/usr/bin/env python3
"""
Calculate quality metrics for 3D chromatin structure.

Metrics include:
- Spearman correlation between distance matrix and contact matrix (dSCC)
- Structure compactness (radius of gyration)
- Local density variations
- Reconstruction loss
"""

import argparse
import json
import sys

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, pearsonr


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


def load_contact_matrix(matrix_file):
    """Load contact matrix from file."""
    return np.loadtxt(matrix_file)


def calculate_dscc(coords, contact_matrix):
    """
    Calculate distance-based Spearman Correlation Coefficient.
    
    This measures how well the 3D distances correlate with contact frequencies.
    Higher contact frequency should correspond to shorter 3D distance.
    """
    if len(coords) == 0:
        return 0.0, 0.0
    
    n = len(coords)
    
    # Truncate to same size
    if contact_matrix.shape[0] > n:
        contact_matrix = contact_matrix[:n, :n]
    elif contact_matrix.shape[0] < n:
        n = contact_matrix.shape[0]
        coords = coords[:n]
    
    # Calculate distance matrix
    dist_matrix = cdist(coords, coords)
    
    # Get upper triangular (excluding diagonal)
    idx = np.triu_indices(n, k=1)
    distances = dist_matrix[idx]
    contacts = contact_matrix[idx]
    
    # Filter out zero contacts
    mask = contacts > 0
    if mask.sum() < 3:
        return 0.0, 0.0
    
    distances_filt = distances[mask]
    contacts_filt = contacts[mask]
    
    # dSCC: correlation between distances and 1/IF (or -IF)
    # Higher IF should mean shorter distance
    spearman_corr, _ = spearmanr(-contacts_filt, distances_filt)
    pearson_corr, _ = pearsonr(-contacts_filt, distances_filt)
    
    return spearman_corr if not np.isnan(spearman_corr) else 0.0, \
           pearson_corr if not np.isnan(pearson_corr) else 0.0


def calculate_radius_of_gyration(coords):
    """Calculate radius of gyration as measure of compactness."""
    if len(coords) == 0:
        return 0.0
    
    center = coords.mean(axis=0)
    rg = np.sqrt(np.mean(np.sum((coords - center)**2, axis=1)))
    return rg


def calculate_local_density(coords, window=5):
    """Calculate local density variations along the structure."""
    if len(coords) < window * 2:
        return 0.0, 0.0
    
    densities = []
    
    for i in range(window, len(coords) - window):
        local_coords = coords[i-window:i+window+1]
        center = local_coords.mean(axis=0)
        local_rg = np.sqrt(np.mean(np.sum((local_coords - center)**2, axis=1)))
        densities.append(1.0 / (local_rg + 1e-6))
    
    return np.mean(densities), np.std(densities)


def calculate_contact_reconstruction_loss(coords, contact_matrix, alpha=1.0):
    """
    Calculate reconstruction loss: how well do distances explain contacts.
    
    Uses the model: contact = 1 / distance^alpha
    """
    if len(coords) == 0:
        return float('inf')
    
    n = len(coords)
    
    if contact_matrix.shape[0] > n:
        contact_matrix = contact_matrix[:n, :n]
    elif contact_matrix.shape[0] < n:
        n = contact_matrix.shape[0]
        coords = coords[:n]
    
    dist_matrix = cdist(coords, coords)
    
    # Predicted contacts from distances
    with np.errstate(divide='ignore', invalid='ignore'):
        predicted = 1.0 / (dist_matrix ** alpha + 1e-6)
    
    # Normalize both matrices
    contact_norm = contact_matrix / (contact_matrix.max() + 1e-6)
    predicted_norm = predicted / (predicted.max() + 1e-6)
    
    # MSE loss
    idx = np.triu_indices(n, k=1)
    loss = np.mean((contact_norm[idx] - predicted_norm[idx])**2)
    
    return loss


def main():
    parser = argparse.ArgumentParser(
        description='Calculate quality metrics for 3D chromatin structure')
    parser.add_argument('--pdb', required=True, help='Input PDB file')
    parser.add_argument('--contact-matrix', required=True,
                       help='Input contact matrix file')
    parser.add_argument('--output', required=True,
                       help='Output JSON file with metrics')
    
    args = parser.parse_args()
    
    metrics = {
        'pdb_file': args.pdb,
        'contact_matrix_file': args.contact_matrix,
        'n_coordinates': 0,
        'dscc_spearman': 0.0,
        'dscc_pearson': 0.0,
        'radius_of_gyration': 0.0,
        'local_density_mean': 0.0,
        'local_density_std': 0.0,
        'reconstruction_loss': float('inf'),
        'status': 'success'
    }
    
    try:
        coords = parse_pdb(args.pdb)
        contact_matrix = load_contact_matrix(args.contact_matrix)
        
        metrics['n_coordinates'] = len(coords)
        
        if len(coords) > 0:
            # dSCC
            dscc_spearman, dscc_pearson = calculate_dscc(coords, contact_matrix)
            metrics['dscc_spearman'] = dscc_spearman
            metrics['dscc_pearson'] = dscc_pearson
            
            # Radius of gyration
            metrics['radius_of_gyration'] = calculate_radius_of_gyration(coords)
            
            # Local density
            ld_mean, ld_std = calculate_local_density(coords)
            metrics['local_density_mean'] = ld_mean
            metrics['local_density_std'] = ld_std
            
            # Reconstruction loss
            metrics['reconstruction_loss'] = calculate_contact_reconstruction_loss(
                coords, contact_matrix)
        else:
            metrics['status'] = 'empty_structure'
            
    except Exception as e:
        metrics['status'] = f'error: {str(e)}'
        print(f"Error: {e}", file=sys.stderr)
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print(f"Metrics written to {args.output}")
    print(f"  dSCC (Spearman): {metrics['dscc_spearman']:.4f}")
    print(f"  Radius of gyration: {metrics['radius_of_gyration']:.2f}")


if __name__ == '__main__':
    main()
