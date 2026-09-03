#!/usr/bin/env python3
"""
Extract small variant statistics from dipcall VCF files.
Outputs per-variant details and aggregated summaries.
"""

import argparse
import gzip
import sys
from collections import defaultdict
import csv


def classify_variant(ref, alt):
    """Classify variant type based on REF and ALT alleles."""
    ref_len = len(ref)
    alt_len = len(alt)
    
    if ref_len == 1 and alt_len == 1:
        return 'SNV'
    elif ref_len > alt_len:
        return 'DEL'
    elif alt_len > ref_len:
        return 'INS'
    elif ref_len == alt_len and ref_len > 1:
        return 'MNV'  # Multi-nucleotide variant
    else:
        return 'COMPLEX'


def get_transition_transversion(ref, alt):
    """Determine if SNV is transition or transversion."""
    if len(ref) != 1 or len(alt) != 1:
        return 'NA'
    
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    ref_upper = ref.upper()
    alt_upper = alt.upper()
    
    if (ref_upper, alt_upper) in transitions:
        return 'transition'
    else:
        return 'transversion'


def parse_vcf_line(line, sample_idx=9):
    """Parse a VCF line and extract relevant variant information."""
    fields = line.strip().split('\t')
    if len(fields) < sample_idx + 1:
        return None
    
    chrom = fields[0]
    pos = int(fields[1])
    var_id = fields[2]
    ref = fields[3]
    alts = fields[4].split(',')
    qual = fields[5]
    filt = fields[6]
    info = fields[7]
    fmt = fields[8]
    sample = fields[sample_idx]
    
    # Parse FORMAT and sample fields
    fmt_fields = fmt.split(':')
    sample_fields = sample.split(':')
    fmt_dict = dict(zip(fmt_fields, sample_fields))
    
    gt = fmt_dict.get('GT', './.')
    
    # Parse allelic depths if available
    ad = fmt_dict.get('AD', '.')
    if ad != '.':
        ad_values = [int(x) if x != '.' else 0 for x in ad.split(',')]
    else:
        ad_values = []
    
    # Determine phasing status
    is_phased = '|' in gt
    
    # Parse genotype alleles
    alleles = gt.replace('|', '/').split('/')
    
    # Handle multi-allelic sites - process each ALT allele
    variants = []
    
    for alt_idx, alt in enumerate(alts, start=1):
        # Check if this allele is present in the genotype
        allele_present = str(alt_idx) in alleles
        
        if not allele_present:
            continue
        
        # Classify variant
        var_type = classify_variant(ref, alt)
        
        # Calculate variant length
        var_len = abs(len(alt) - len(ref))
        
        # Determine zygosity for this allele
        allele_count = alleles.count(str(alt_idx))
        if allele_count == 2:
            zygosity = 'HOM_ALT'
        elif allele_count == 1:
            zygosity = 'HET'
        else:
            zygosity = 'MISSING'
        
        # Determine haplotype
        if is_phased:
            hap1_has = alleles[0] == str(alt_idx) if len(alleles) > 0 else False
            hap2_has = alleles[1] == str(alt_idx) if len(alleles) > 1 else False
            if hap1_has and hap2_has:
                haplotype = 'both'
            elif hap1_has:
                haplotype = 'haplotype1'
            elif hap2_has:
                haplotype = 'haplotype2'
            else:
                haplotype = 'none'
        else:
            haplotype = 'unphased'
        
        # Ti/Tv for SNVs
        ti_tv = get_transition_transversion(ref, alt) if var_type == 'SNV' else 'NA'
        
        # Biallelic check
        is_biallelic = len(alts) == 1
        
        # Allele frequency (from this sample's genotype)
        total_alleles = len([a for a in alleles if a != '.'])
        alt_alleles = alleles.count(str(alt_idx))
        allele_freq = alt_alleles / total_alleles if total_alleles > 0 else 0
        
        variants.append({
            'chrom': chrom,
            'pos': pos,
            'var_id': var_id if var_id != '.' else f"{chrom}_{pos}",
            'ref': ref,
            'alt': alt,
            'var_type': var_type,
            'var_len': var_len,
            'genotype': gt,
            'is_phased': is_phased,
            'zygosity': zygosity,
            'haplotype': haplotype,
            'filter': filt,
            'ti_tv': ti_tv,
            'is_biallelic': is_biallelic,
            'allele_freq': allele_freq
        })
    
    return variants


def process_vcf(vcf_file, sample_name):
    """Process VCF file and extract variant statistics."""
    all_variants = []
    
    opener = gzip.open if vcf_file.endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variants = parse_vcf_line(line)
            if variants:
                for v in variants:
                    v['sample'] = sample_name
                    all_variants.append(v)
    
    return all_variants


def calculate_summary(variants, sample_name):
    """Calculate summary statistics from variants."""
    summary = defaultdict(lambda: defaultdict(int))
    length_sums = defaultdict(lambda: defaultdict(int))
    
    # Count transitions and transversions for Ti/Tv ratio
    ti_count = 0
    tv_count = 0
    
    # Per-type breakdowns
    per_type_zygosity = defaultdict(lambda: defaultdict(int))
    per_type_haplotype = defaultdict(lambda: defaultdict(int))
    per_type_phasing = defaultdict(lambda: defaultdict(int))
    per_type_allelic = defaultdict(lambda: defaultdict(int))
    per_type_indel_size = defaultdict(lambda: defaultdict(int))
    per_type_chrom = defaultdict(lambda: defaultdict(int))
    type_totals = defaultdict(int)
    
    for v in variants:
        var_type = v['var_type']
        type_totals[var_type] += 1
        
        # Count by type
        summary['by_type'][var_type] += 1
        length_sums['by_type'][var_type] += v['var_len']
        
        # Count by zygosity (overall and per-type)
        summary['by_zygosity'][v['zygosity']] += 1
        per_type_zygosity[var_type][v['zygosity']] += 1
        
        # Count by haplotype (overall and per-type)
        summary['by_haplotype'][v['haplotype']] += 1
        per_type_haplotype[var_type][v['haplotype']] += 1
        
        # Count by phasing status (overall and per-type)
        phasing_key = 'phased' if v['is_phased'] else 'unphased'
        summary['by_phasing'][phasing_key] += 1
        per_type_phasing[var_type][phasing_key] += 1
        
        # Count by filter status
        summary['by_filter'][v['filter']] += 1
        
        # Count biallelic vs multiallelic (overall and per-type)
        biallelic_key = 'biallelic' if v['is_biallelic'] else 'multiallelic'
        summary['by_allelic'][biallelic_key] += 1
        per_type_allelic[var_type][biallelic_key] += 1
        
        # Count Ti/Tv
        if v['ti_tv'] == 'transition':
            ti_count += 1
        elif v['ti_tv'] == 'transversion':
            tv_count += 1
        
        # Count by chromosome (overall and per-type for major types)
        summary['by_chromosome'][v['chrom']] += 1
        if var_type in ['SNV', 'DEL', 'INS']:
            per_type_chrom[var_type][v['chrom']] += 1
        
        # Size distribution for indels (overall and per-type)
        if var_type in ['DEL', 'INS']:
            if v['var_len'] == 1:
                size_class = '1bp'
            elif v['var_len'] <= 5:
                size_class = '2-5bp'
            elif v['var_len'] <= 10:
                size_class = '6-10bp'
            elif v['var_len'] <= 20:
                size_class = '11-20bp'
            elif v['var_len'] <= 50:
                size_class = '21-50bp'
            else:
                size_class = '>50bp'
            summary['indel_size'][size_class] += 1
            per_type_indel_size[var_type][size_class] += 1
    
    # Create summary rows
    summary_rows = []
    
    # Total count
    total_variants = len(variants)
    summary_rows.append({
        'sample': sample_name,
        'category': 'total',
        'subcategory': 'all_variants',
        'var_type': 'all',
        'count': total_variants,
        'percentage': 100.0
    })
    
    # By type with percentages
    for var_type, count in sorted(summary['by_type'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'variant_type',
            'subcategory': var_type,
            'var_type': var_type,
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # Ti/Tv ratio
    ti_tv_ratio = ti_count / tv_count if tv_count > 0 else 0
    summary_rows.append({
        'sample': sample_name,
        'category': 'ti_tv',
        'subcategory': 'transitions',
        'var_type': 'SNV',
        'count': ti_count,
        'percentage': ti_tv_ratio  # Store ratio in percentage field
    })
    summary_rows.append({
        'sample': sample_name,
        'category': 'ti_tv',
        'subcategory': 'transversions',
        'var_type': 'SNV',
        'count': tv_count,
        'percentage': ti_tv_ratio
    })
    
    # By zygosity (overall)
    for zyg, count in sorted(summary['by_zygosity'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'zygosity',
            'subcategory': zyg,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # By zygosity (per variant type)
    for var_type in sorted(per_type_zygosity.keys()):
        type_total = type_totals[var_type]
        for zyg, count in sorted(per_type_zygosity[var_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'zygosity',
                'subcategory': zyg,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    # By haplotype (overall)
    for hap, count in sorted(summary['by_haplotype'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'haplotype',
            'subcategory': hap,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # By haplotype (per variant type)
    for var_type in sorted(per_type_haplotype.keys()):
        type_total = type_totals[var_type]
        for hap, count in sorted(per_type_haplotype[var_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'haplotype',
                'subcategory': hap,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    # By phasing (overall)
    for phasing, count in sorted(summary['by_phasing'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'phasing',
            'subcategory': phasing,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # By phasing (per variant type)
    for var_type in sorted(per_type_phasing.keys()):
        type_total = type_totals[var_type]
        for phasing, count in sorted(per_type_phasing[var_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'phasing',
                'subcategory': phasing,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    # By allelic status (overall)
    for allelic, count in sorted(summary['by_allelic'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'allelic',
            'subcategory': allelic,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # By allelic status (per variant type)
    for var_type in sorted(per_type_allelic.keys()):
        type_total = type_totals[var_type]
        for allelic, count in sorted(per_type_allelic[var_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'allelic',
                'subcategory': allelic,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    # By filter (overall only)
    for filt, count in sorted(summary['by_filter'].items()):
        filt_name = filt if filt != '.' else 'PASS'
        summary_rows.append({
            'sample': sample_name,
            'category': 'filter',
            'subcategory': filt_name,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # Indel size distribution (overall)
    for size_class, count in sorted(summary['indel_size'].items()):
        summary_rows.append({
            'sample': sample_name,
            'category': 'indel_size',
            'subcategory': size_class,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # Indel size distribution (per variant type: DEL and INS)
    for var_type in sorted(per_type_indel_size.keys()):
        type_total = type_totals[var_type]
        for size_class, count in sorted(per_type_indel_size[var_type].items()):
            summary_rows.append({
                'sample': sample_name,
                'category': 'indel_size',
                'subcategory': size_class,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    # By chromosome (overall)
    for chrom, count in sorted(summary['by_chromosome'].items(), 
                                key=lambda x: (x[0].replace('chr', '').zfill(2))):
        summary_rows.append({
            'sample': sample_name,
            'category': 'chromosome',
            'subcategory': chrom,
            'var_type': 'all',
            'count': count,
            'percentage': (count / total_variants * 100) if total_variants else 0
        })
    
    # By chromosome (per variant type: SNV, DEL, INS)
    for var_type in sorted(per_type_chrom.keys()):
        type_total = type_totals[var_type]
        for chrom, count in sorted(per_type_chrom[var_type].items(),
                                    key=lambda x: (x[0].replace('chr', '').zfill(2))):
            summary_rows.append({
                'sample': sample_name,
                'category': 'chromosome',
                'subcategory': chrom,
                'var_type': var_type,
                'count': count,
                'percentage': (count / type_total * 100) if type_total else 0
            })
    
    return summary_rows


def main():
    parser = argparse.ArgumentParser(description='Extract small variant statistics from dipcall VCF')
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
    detail_fields = ['sample', 'chrom', 'pos', 'var_id', 'ref', 'alt', 'var_type', 
                     'var_len', 'genotype', 'is_phased', 'zygosity', 'haplotype', 
                     'filter', 'ti_tv', 'is_biallelic', 'allele_freq']
    
    with open(args.output_details, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=detail_fields, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(variants)
    
    # Calculate and write summary
    summary = calculate_summary(variants, args.sample)
    
    summary_fields = ['sample', 'category', 'subcategory', 'var_type', 'count', 'percentage']
    
    with open(args.output_summary, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
    
    print(f"Processed {len(variants)} variants from {args.sample}", file=sys.stderr)


if __name__ == '__main__':
    main()
