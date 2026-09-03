#!/usr/bin/env python3
"""
Extract structural variant statistics from hapdiff VCF files.
Outputs per-variant details and aggregated summaries.
"""

import argparse
import gzip
import sys
from collections import defaultdict
import csv


def parse_vcf_line(line):
    """Parse a VCF line and extract relevant SV information."""
    fields = line.strip().split('\t')
    if len(fields) < 10:
        return None
    
    chrom = fields[0]
    pos = int(fields[1])
    sv_id = fields[2]
    ref = fields[3]
    alt = fields[4]
    qual = fields[5]
    filt = fields[6]
    info = fields[7]
    fmt = fields[8]
    sample = fields[9]
    
    # Parse INFO field
    info_dict = {}
    for item in info.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    # Extract SV type
    sv_type = info_dict.get('SVTYPE', 'UNKNOWN')
    
    # Extract SV length
    sv_len = abs(int(info_dict.get('SVLEN', 0)))
    
    # If SVLEN is 0, calculate from REF/ALT
    if sv_len == 0:
        if sv_type == 'DEL':
            sv_len = len(ref) - len(alt)
        elif sv_type == 'INS':
            sv_len = len(alt) - len(ref)
    
    # Extract END position
    end = int(info_dict.get('END', pos))
    
    # Parse genotype
    fmt_fields = fmt.split(':')
    sample_fields = sample.split(':')
    fmt_dict = dict(zip(fmt_fields, sample_fields))
    
    gt = fmt_dict.get('GT', './.')
    
    # Determine phasing status
    is_phased = '|' in gt
    
    # Determine zygosity
    alleles = gt.replace('|', '/').split('/')
    if alleles[0] == alleles[1]:
        if alleles[0] == '0':
            zygosity = 'REF_HOM'
        else:
            zygosity = 'ALT_HOM'
    else:
        zygosity = 'HET'
    
    # Determine which haplotype has the variant
    if is_phased:
        hap1_alt = alleles[0] != '0' and alleles[0] != '.'
        hap2_alt = alleles[1] != '0' and alleles[1] != '.'
        if hap1_alt and hap2_alt:
            haplotype = 'both'
        elif hap1_alt:
            haplotype = 'haplotype1'
        elif hap2_alt:
            haplotype = 'haplotype2'
        else:
            haplotype = 'none'
    else:
        haplotype = 'unphased'
    
    # Classify SV by size
    if sv_len < 50:
        size_class = 'small (<50bp)'
    elif sv_len < 500:
        size_class = 'medium (50-500bp)'
    elif sv_len < 5000:
        size_class = 'large (500bp-5kb)'
    elif sv_len < 50000:
        size_class = 'very_large (5kb-50kb)'
    else:
        size_class = 'massive (>50kb)'
    
    return {
        'chrom': chrom,
        'pos': pos,
        'end': end,
        'sv_id': sv_id,
        'sv_type': sv_type,
        'sv_len': sv_len,
        'size_class': size_class,
        'genotype': gt,
        'is_phased': is_phased,
        'zygosity': zygosity,
        'haplotype': haplotype,
        'filter': filt
    }


def process_vcf(vcf_file, sample_name):
    """Process VCF file and extract SV statistics."""
    variants = []
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            sv_info = parse_vcf_line(line)
            if sv_info:
                sv_info['sample'] = sample_name
                variants.append(sv_info)
    
    return variants


def calculate_summary(variants, sample_name):
    """Calculate summary statistics from variants."""
    summary = defaultdict(lambda: defaultdict(int))
    length_sums = defaultdict(lambda: defaultdict(int))
    
    # Per-type breakdowns
    per_type_zygosity = defaultdict(lambda: defaultdict(int))
    per_type_haplotype = defaultdict(lambda: defaultdict(int))
    per_type_size_class = defaultdict(lambda: defaultdict(int))
    per_type_chrom = defaultdict(lambda: defaultdict(int))
    
    for v in variants:
        sv_type = v['sv_type']
        
        # Count by type
        summary['by_type'][sv_type] += 1
        length_sums['by_type'][sv_type] += v['sv_len']
        
        # Count by zygosity (overall and per-type)
        summary['by_zygosity'][v['zygosity']] += 1
        per_type_zygosity[sv_type][v['zygosity']] += 1
        
        # Count by size class (overall and per-type)
        summary['by_size_class'][v['size_class']] += 1
        per_type_size_class[sv_type][v['size_class']] += 1
        
        # Count by haplotype (overall and per-type)
        summary['by_haplotype'][v['haplotype']] += 1
        per_type_haplotype[sv_type][v['haplotype']] += 1
        
        # Count by phasing status
        phasing_key = 'phased' if v['is_phased'] else 'unphased'
        summary['by_phasing'][phasing_key] += 1
        
        # Count by filter status
        summary['by_filter'][v['filter']] += 1
        
        # Count by chromosome (overall and per-type)
        summary['by_chromosome'][v['chrom']] += 1
        per_type_chrom[sv_type][v['chrom']] += 1
    
    # Create summary rows
    summary_rows = []
    
    # Total count
    summary_rows.append({
        'sample': sample_name,
        'category': 'total',
        'subcategory': 'all_svs',
        'sv_type': 'all',
        'count': len(variants),
        'total_length': sum(v['sv_len'] for v in variants),
        'mean_length': sum(v['sv_len'] for v in variants) / len(variants) if variants else 0
    })
    
    # By type
    for sv_type, count in sorted(summary['by_type'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'sv_type',
            'subcategory': sv_type,
            'sv_type': sv_type,
            'count': count,
            'total_length': length_sums['by_type'][sv_type],
            'mean_length': length_sums['by_type'][sv_type] / count if count else 0
        })
    
    # By zygosity (overall)
    for zyg, count in sorted(summary['by_zygosity'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'zygosity',
            'subcategory': zyg,
            'sv_type': 'all',
            'count': count,
            'total_length': 0,
            'mean_length': 0
        })
    
    # By zygosity (per SV type)
    for sv_type in sorted(per_type_zygosity.keys()):
        for zyg, count in sorted(per_type_zygosity[sv_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'zygosity',
                'subcategory': zyg,
                'sv_type': sv_type,
                'count': count,
                'total_length': 0,
                'mean_length': 0
            })
    
    # By size class (overall)
    for size_class, count in sorted(summary['by_size_class'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'size_class',
            'subcategory': size_class,
            'sv_type': 'all',
            'count': count,
            'total_length': 0,
            'mean_length': 0
        })
    
    # By size class (per SV type)
    for sv_type in sorted(per_type_size_class.keys()):
        for size_class, count in sorted(per_type_size_class[sv_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'size_class',
                'subcategory': size_class,
                'sv_type': sv_type,
                'count': count,
                'total_length': 0,
                'mean_length': 0
            })
    
    # By haplotype (overall)
    for hap, count in sorted(summary['by_haplotype'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'haplotype',
            'subcategory': hap,
            'sv_type': 'all',
            'count': count,
            'total_length': 0,
            'mean_length': 0
        })
    
    # By haplotype (per SV type)
    for sv_type in sorted(per_type_haplotype.keys()):
        for hap, count in sorted(per_type_haplotype[sv_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'haplotype',
                'subcategory': hap,
                'sv_type': sv_type,
                'count': count,
                'total_length': 0,
                'mean_length': 0
            })
    
    # By phasing (overall only - same for all types)
    for phasing, count in sorted(summary['by_phasing'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'phasing',
            'subcategory': phasing,
            'sv_type': 'all',
            'count': count,
            'total_length': 0,
            'mean_length': 0
        })
    
    # By chromosome (overall)
    for chrom, count in sorted(summary['by_chromosome'].items(), 
                                key=lambda x: (x[0].replace('chr', '').zfill(2))):
        summary_rows.append({
            'sample': sample_name,
            'category': 'chromosome',
            'subcategory': chrom,
            'sv_type': 'all',
            'count': count,
            'total_length': 0,
            'mean_length': 0
        })
    
    # By chromosome (per SV type - only for major types DEL, INS)
    for sv_type in ['DEL', 'INS']:
        if sv_type in per_type_chrom:
            for chrom, count in sorted(per_type_chrom[sv_type].items(), 
                                        key=lambda x: (x[0].replace('chr', '').zfill(2))):
                summary_rows.append({
                    'sample': sample_name,
                    'category': 'chromosome',
                    'subcategory': chrom,
                    'sv_type': sv_type,
                    'count': count,
                    'total_length': 0,
                    'mean_length': 0
                })
    
    return summary_rows


def main():
    parser = argparse.ArgumentParser(description='Extract SV statistics from hapdiff VCF')
    parser.add_argument('--vcf', required=True, help='Input VCF file (gzipped supported)')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--output-details', required=True, help='Output TSV with per-variant details')
    parser.add_argument('--output-summary', required=True, help='Output TSV with summary statistics')
    
    args = parser.parse_args()
    
    # Process VCF
    variants = process_vcf(args.vcf, args.sample)
    
    if not variants:
        print(f"Warning: No variants found in {args.vcf}", file=sys.stderr)
    
    # Write detailed output
    detail_fields = ['sample', 'chrom', 'pos', 'end', 'sv_id', 'sv_type', 'sv_len', 
                     'size_class', 'genotype', 'is_phased', 'zygosity', 'haplotype', 'filter']
    
    with open(args.output_details, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=detail_fields, delimiter='\t', 
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(variants)
    
    # Calculate and write summary
    summary = calculate_summary(variants, args.sample)
    
    summary_fields = ['sample', 'category', 'subcategory', 'sv_type', 'count', 'total_length', 'mean_length']
    
    with open(args.output_summary, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
    
    print(f"Processed {len(variants)} SVs from {args.sample}", file=sys.stderr)


if __name__ == '__main__':
    main()
