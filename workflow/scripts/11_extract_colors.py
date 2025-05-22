#!/usr/bin/env python
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Extract colors from PAF file')
    parser.add_argument('-i', '--input', required=True, help='Input PAF file')
    parser.add_argument('-c', '--colours', required=False, help='Input colours from Phasing (Optional)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-l', '--legend', action='store_true', help='Generate color legend')
    return parser.parse_args()


chromosome_colors = {
    # Main chromosomes - using colorblind-friendly base colors with light/dark variants
    'chr1':  {'haplotype1': '#E69F00', 'haplotype2': '#FFB319'},
    'chr2':  {'haplotype1': '#56B4E9', 'haplotype2': '#7BC6EE'},
    'chr3':  {'haplotype1': '#009E73', 'haplotype2': '#00BF8C'},
    'chr4':  {'haplotype1': '#CC79A7', 'haplotype2': '#D694B8'},
    'chr5':  {'haplotype1': '#0072B2', 'haplotype2': '#0089D9'},
    'chr6':  {'haplotype1': '#D55E00', 'haplotype2': '#FF7400'},
    'chr7':  {'haplotype1': '#666666', 'haplotype2': '#999999'},
    
    # Recycling colors with different shades
    'chr8':  {'haplotype1': '#E6AB23', 'haplotype2': '#FFCD66'},
    'chr9':  {'haplotype1': '#5699E9', 'haplotype2': '#89B9F0'},
    'chr10': {'haplotype1': '#00A880', 'haplotype2': '#33BF99'},
    'chr11': {'haplotype1': '#CC8DB3', 'haplotype2': '#D6A8C4'},
    'chr12': {'haplotype1': '#1A7FBA', 'haplotype2': '#4D9ECC'},
    'chr13': {'haplotype1': '#D57533', 'haplotype2': '#FF9B66'},
    'chr14': {'haplotype1': '#737373', 'haplotype2': '#A6A6A6'},
    'chr15': {'haplotype1': '#E6B847', 'haplotype2': '#FFD480'},
    'chr16': {'haplotype1': '#567DE9', 'haplotype2': '#89A3F0'},
    'chr17': {'haplotype1': '#00B28C', 'haplotype2': '#33C6A6'},
    'chr18': {'haplotype1': '#CCA0BF', 'haplotype2': '#D6BBD0'},
    'chr19': {'haplotype1': '#338CC2', 'haplotype2': '#66ACD9'},
    'chr20': {'haplotype1': '#D58C66', 'haplotype2': '#FFB499'},
    'chr21': {'haplotype1': '#808080', 'haplotype2': '#B3B3B3'},
    'chr22': {'haplotype1': '#E6C47F', 'haplotype2': '#FFD699'},

    # Sex chromosomes
    'chrX':  {'haplotype1': '#9467BD', 'haplotype2': '#B189D6'},
    'chrY':  {'haplotype1': '#8C564B', 'haplotype2': '#A67C73'},

     # Mitochondrial chromosome
    'chrM': {'haplotype1': '#17BECF', 'haplotype2': '#4DD8E6'}
}

def plot_color_legend(chromosome_colors, output_file):
    """
    Generate a color legend visualization for chromosome colors in a single column.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches  
    import numpy as np
    
    # Sort chromosomes for better display
    sorted_chromosomes = sorted(chromosome_colors.keys(), 
                               key=lambda x: int(x[3:]) if x[3:].isdigit() else 
                                            (100 if x == 'chrX' else 
                                             101 if x == 'chrY' else 
                                             102 if x == 'chrM' else 999))
    
    # Set up the figure with appropriate size for a single column
    num_chromosomes = len(sorted_chromosomes)
    fig_height = max(10, num_chromosomes * 0.4)  # Adjust height based on number of chromosomes
    fig, ax = plt.subplots(figsize=(8, fig_height))
    
    # Remove axes
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    
    # Calculate spacing
    row_height = 0.9 / num_chromosomes
    start_y = 0.95
    
    # Add legend items
    for i, chrom in enumerate(sorted_chromosomes):
        y_pos = start_y - i * row_height
        
        # Add chromosome name
        ax.text(0.05, y_pos, chrom, fontsize=12, fontweight='bold')
        
        # Add color patches for haplotypes
        hap1_color = chromosome_colors[chrom]['haplotype1']
        hap2_color = chromosome_colors[chrom]['haplotype2']
        
        # Haplotype 1 patch - smaller bar
        bar_height = row_height * 0.4
        bar_width = 0.2
        
        # Haplotype 1 patch
        hap1_patch = mpatches.Rectangle((0.3, y_pos - bar_height/2), bar_width, bar_height, 
                                        facecolor=hap1_color, edgecolor='black')
        ax.add_patch(hap1_patch)
        ax.text(0.3 + bar_width + 0.02, y_pos, 'Hap 1', fontsize=10, va='center')
        
        # Haplotype 2 patch
        hap2_patch = mpatches.Rectangle((0.65, y_pos - bar_height/2), bar_width, bar_height, 
                                        facecolor=hap2_color, edgecolor='black')
        ax.add_patch(hap2_patch)
        ax.text(0.65 + bar_width + 0.02, y_pos, 'Hap 2', fontsize=10, va='center')
    
    # Add title
    plt.suptitle('Chromosome Color Legend', fontsize=16, y=0.98)
    
    # Add subtitle explaining the color scheme
    plt.figtext(0.5, 0.01, 
                'Color scheme: Colorblind-friendly palette with light/dark variants for haplotypes', 
                ha='center', fontsize=10)
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Color legend saved to {output_file}")

def main():
    args = parse_args()
    
    # Dictionary to store best match for each contig
    contig_best_matches = {}
    
    with open(args.input, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            contig_id = fields[0]
            chrom = fields[5]
            match_q = int(fields[10])  # Match length
            match_l = int(fields[11])  # Match quality
            

            # If this is the first time we see this contig or if this match is better
            if contig_id not in contig_best_matches or \
               (match_q >= contig_best_matches[contig_id]['match_q'] and 
                match_l >= contig_best_matches[contig_id]['match_l']):
                contig_best_matches[contig_id] = {
                    'chrom': chrom,
                    'match_q': match_q,
                    'match_l': match_l
                }
    
    # Generate the output lines
    contig_colors = []
    for contig_id, match_info in contig_best_matches.items():
        chrom = match_info['chrom']
        color = chromosome_colors[chrom]['haplotype1']
        contig_colors.append(f"{contig_id}\t{color}\t{chrom}")
    
    with open(args.output, "w") as f:
        f.write("contig\tcolor\tchromosome\n")
        f.write("\n".join(contig_colors))

    if args.legend:
        plot_color_legend(chromosome_colors, os.path.join(os.path.dirname(args.output), 'chromosome_color_legend.png'))

if __name__ == "__main__":
    main()