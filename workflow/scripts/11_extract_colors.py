#!/usr/bin/env python
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Extract colors from PAF file')
    parser.add_argument('-i', '--input', required=True, help='Input PAF file')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    parser.add_argument('-p', '--haplotype', default="haplotype1", 
                        choices=["haplotype1", "haplotype2", "unphased"],
                        help='Haplotype to use for coloring (default: haplotype1)')
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

def main():
    args = parse_args()
    
    # Dictionary to store best match for each contig
    contig_best_matches = {}
    
    # Set haplotype - treat "unphased" as "haplotype1"
    hp = args.haplotype
    if hp == "unphased":
        hp = "haplotype1"
    
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
        color = chromosome_colors[chrom][hp]
        contig_colors.append(f"{contig_id}\t{color}")
    
    with open(args.output, "w") as f:
        f.write("contig\tcolor\n")
        f.write("\n".join(contig_colors))

if __name__ == "__main__":
    main()