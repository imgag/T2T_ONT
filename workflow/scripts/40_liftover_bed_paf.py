#!/usr/bin/env python3
"""
Liftover BED files using PAF alignments while preserving all original columns.
"""

import argparse
import sys
import re
from collections import defaultdict

class PAFAlignment:
    def __init__(self, line):
        fields = line.strip().split('\t')
        self.query_name = fields[0]
        self.query_len = int(fields[1])
        self.query_start = int(fields[2])
        self.query_end = int(fields[3])
        self.strand = fields[4]
        self.target_name = fields[5]
        self.target_len = int(fields[6])
        self.target_start = int(fields[7])
        self.target_end = int(fields[8])
        self.matches = int(fields[9])
        self.alignment_len = int(fields[10])
        self.quality = int(fields[11])
        
        # Parse CIGAR string from cg tag
        self.cigar = None
        for field in fields[12:]:
            if field.startswith('cg:Z:'):
                self.cigar = field[5:]
                break

class BEDEntry:
    def __init__(self, line):
        fields = line.strip().split('\t')
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.name = fields[3] if len(fields) > 3 else "."
        self.score = fields[4] if len(fields) > 4 else "0"
        self.strand = fields[5] if len(fields) > 5 else "+"
        # Preserve additional BED columns (thickStart, thickEnd, itemRgb, etc.)
        self.extra_fields = fields[6:] if len(fields) > 6 else []
        self.original_line = line.strip()

def load_paf_alignments(paf_file, min_mapq=5, min_len=50000):
    """Load PAF alignments and create coordinate mapping."""
    alignments = defaultdict(list)
    filtered_stats = {
        'total': 0,
        'low_mapq': 0,
        'short_len': 0,
        'passed': 0
    }
    
    with open(paf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            filtered_stats['total'] += 1
            aln = PAFAlignment(line)
            
            # Track filtering reasons
            if aln.quality < min_mapq:
                filtered_stats['low_mapq'] += 1
                continue
            if aln.alignment_len < min_len:
                filtered_stats['short_len'] += 1
                continue
            
            filtered_stats['passed'] += 1
            alignments[aln.query_name].append(aln)
    
    # Sort alignments by query coordinates
    for query_name in alignments:
        alignments[query_name].sort(key=lambda x: x.query_start)
    
    print(f"PAF alignment filtering stats:", file=sys.stderr)
    print(f"  Total alignments: {filtered_stats['total']}", file=sys.stderr)
    print(f"  Filtered (low MAPQ < {min_mapq}): {filtered_stats['low_mapq']}", file=sys.stderr)
    print(f"  Filtered (short length < {min_len}): {filtered_stats['short_len']}", file=sys.stderr)
    print(f"  Passed filters: {filtered_stats['passed']}", file=sys.stderr)
    
    return alignments

def parse_cigar_operations(cigar):
    """Parse CIGAR string into operations."""
    operations = []
    for match in re.finditer(r'(\d+)([MID])', cigar):
        length = int(match.group(1))
        op = match.group(2)
        operations.append((length, op))
    return operations

def liftover_position(pos, aln, debug=False):
    """
    Liftover a single position using CIGAR operations.
    Returns target coordinate or None if position falls in deletion.
    """
    if not aln.cigar:
        # Simple linear mapping without CIGAR
        if aln.strand == '+':
            offset = pos - aln.query_start
            result = aln.target_start + offset
        else:
            offset = aln.query_end - pos
            result = aln.target_start + offset
        
        if debug:
            print(f"    Linear mapping: pos={pos}, offset={offset}, result={result}", file=sys.stderr)
        return result
    
    # Use CIGAR for precise mapping
    operations = parse_cigar_operations(aln.cigar)
    
    query_pos = aln.query_start
    target_pos = aln.target_start
    
    if debug:
        print(f"    CIGAR mapping: pos={pos}, operations={operations[:5]}...", file=sys.stderr)
    
    for length, op in operations:
        if op == 'M':  # Match/mismatch
            if query_pos <= pos < query_pos + length:
                if aln.strand == '+':
                    result = target_pos + (pos - query_pos)
                else:
                    result = target_pos + length - (pos - query_pos) - 1
                if debug:
                    print(f"    Found in match block: query_pos={query_pos}, target_pos={target_pos}, result={result}", file=sys.stderr)
                return result
            query_pos += length
            target_pos += length
        elif op == 'I':  # Insertion in query
            if query_pos <= pos < query_pos + length:
                if debug:
                    print(f"    Position falls in insertion (query {query_pos}-{query_pos+length})", file=sys.stderr)
                return None  # Position falls in insertion (not in target)
            query_pos += length
        elif op == 'D':  # Deletion in query
            target_pos += length
    
    if debug:
        print(f"    Position {pos} not found in any CIGAR block", file=sys.stderr)
    return None

def liftover_bed_region(bed_entry, alignments, debug=False):
    """
    Liftover a BED region, returning list of lifted regions.
    """
    failure_reasons = []
    
    if bed_entry.chrom not in alignments:
        reason = f"No alignments found for chromosome '{bed_entry.chrom}'"
        failure_reasons.append(reason)
        if debug:
            print(f"  FAIL: {reason}", file=sys.stderr)
            print(f"        Available chromosomes: {list(alignments.keys())[:10]}...", file=sys.stderr)
        return [], failure_reasons
    
    lifted_regions = []
    overlap_count = 0
    
    for i, aln in enumerate(alignments[bed_entry.chrom]):
        # Check if region overlaps with alignment
        overlap_start = max(bed_entry.start, aln.query_start)
        overlap_end = min(bed_entry.end, aln.query_end)
        
        if overlap_start >= overlap_end:
            if debug and i < 3:  # Only show first few for brevity
                print(f"    Alignment {i}: no overlap", file=sys.stderr)
                print(f"      BED region: {bed_entry.start}-{bed_entry.end}", file=sys.stderr)
                print(f"      Alignment: {aln.query_start}-{aln.query_end} (MAPQ={aln.quality})", file=sys.stderr)
            continue  # No overlap
        
        overlap_count += 1
        if debug:
            print(f"    Alignment {i}: overlap found", file=sys.stderr)
            print(f"      BED region: {bed_entry.start}-{bed_entry.end}", file=sys.stderr)
            print(f"      Alignment: {aln.query_start}-{aln.query_end} -> {aln.target_start}-{aln.target_end}", file=sys.stderr)
            print(f"      Overlap: {overlap_start}-{overlap_end}", file=sys.stderr)
        
        # Liftover start and end positions
        lifted_start = liftover_position(overlap_start, aln, debug)
        lifted_end = liftover_position(overlap_end - 1, aln, debug)  # BED end is exclusive
        
        if lifted_start is None:
            reason = f"Start position {overlap_start} falls in insertion/unmappable region"
            failure_reasons.append(reason)
            if debug:
                print(f"      FAIL: {reason}", file=sys.stderr)
            continue
        
        if lifted_end is None:
            reason = f"End position {overlap_end-1} falls in insertion/unmappable region"
            failure_reasons.append(reason)
            if debug:
                print(f"      FAIL: {reason}", file=sys.stderr)
            continue
        
        # Ensure proper ordering
        if lifted_start > lifted_end:
            lifted_start, lifted_end = lifted_end, lifted_start
        
        lifted_end += 1  # Convert back to exclusive end
        
        if debug:
            print(f"      SUCCESS: {overlap_start}-{overlap_end} -> {lifted_start}-{lifted_end}", file=sys.stderr)
        
        # Create lifted BED entry preserving all original fields
        lifted_name = f"{bed_entry.chrom}_{overlap_start}_{overlap_end}"
        if overlap_start > bed_entry.start or overlap_end < bed_entry.end:
            lifted_name += "_partial"
        
        # Preserve original name and add liftover info
        if bed_entry.name != ".":
            lifted_name = bed_entry.name
        
        lifted_regions.append({
            'chrom': aln.target_name,
            'start': lifted_start,
            'end': lifted_end,
            'name': lifted_name,
            'score': bed_entry.score,
            'strand': aln.strand,
            'extra_fields': bed_entry.extra_fields,
            'original_entry': bed_entry
        })
    
    if not lifted_regions and overlap_count == 0:
        reason = f"No overlapping alignments found for region {bed_entry.start}-{bed_entry.end}"
        failure_reasons.append(reason)
        if debug:
            # Show nearby alignments for context
            nearby_alns = [aln for aln in alignments[bed_entry.chrom] 
                          if abs(aln.query_start - bed_entry.start) < 1000000 or 
                             abs(aln.query_end - bed_entry.end) < 1000000][:3]
            print(f"  FAIL: {reason}", file=sys.stderr)
            print(f"        Nearby alignments:", file=sys.stderr)
            for aln in nearby_alns:
                print(f"          {aln.query_start}-{aln.query_end} (MAPQ={aln.quality})", file=sys.stderr)
    
    return lifted_regions, failure_reasons

def liftover_bed_file(bed_file, paf_file, output_file, unmapped_file, debug=False):
    """Main liftover function."""
    
    print(f"Loading PAF alignments from {paf_file}", file=sys.stderr)
    alignments = load_paf_alignments(paf_file)
    
    print(f"Found alignments for {len(alignments)} query sequences", file=sys.stderr)
    
    lifted_count = 0
    unmapped_count = 0
    failure_summary = defaultdict(int)
    
    with open(bed_file, 'r') as infile, \
         open(output_file, 'w') as outfile, \
         open(unmapped_file, 'w') as unmapped:
        
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if line.startswith('#') or line.startswith('track') or not line:
                # Preserve header lines in output
                if line.startswith('track'):
                    outfile.write(line + '\n')
                continue
            
            bed_entry = BEDEntry(line)
            
            if debug:
                print(f"\nProcessing BED entry {line_num}: {bed_entry.chrom}:{bed_entry.start}-{bed_entry.end} ({bed_entry.name})", file=sys.stderr)
            
            lifted_regions, failure_reasons = liftover_bed_region(bed_entry, alignments, debug)
            
            if not lifted_regions:
                unmapped.write(line + '\n')
                
                # Write detailed failure reason to unmapped file
                if failure_reasons:
                    unmapped.write(f"# FAILURE REASONS: {'; '.join(failure_reasons)}\n")
                    for reason in failure_reasons:
                        failure_summary[reason] += 1
                
                unmapped_count += 1
                continue
            
            # Write lifted regions, preserving all original BED columns
            for region in lifted_regions:
                output_fields = [
                    region['chrom'],
                    str(region['start']),
                    str(region['end']),
                    region['name'],
                    region['score'],
                    region['strand']
                ]
                
                # Add extra fields (thickStart, thickEnd, itemRgb, etc.)
                output_fields.extend(region['extra_fields'])
                
                outfile.write('\t'.join(output_fields) + '\n')
                lifted_count += 1
    
    print(f"\nSUMMARY:", file=sys.stderr)
    print(f"Successfully lifted {lifted_count} regions", file=sys.stderr)
    print(f"Failed to lift {unmapped_count} regions", file=sys.stderr)
    
    if failure_summary:
        print(f"\nFailure reasons:", file=sys.stderr)
        for reason, count in sorted(failure_summary.items(), key=lambda x: x[1], reverse=True):
            print(f"  {count:4d}x: {reason}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Liftover BED files using PAF alignments')
    parser.add_argument('--bed', required=True, help='Input BED file')
    parser.add_argument('--paf', required=True, help='PAF alignment file')
    parser.add_argument('--output', required=True, help='Output lifted BED file')
    parser.add_argument('--unmapped', required=True, help='Output unmapped regions file')
    parser.add_argument('--min-mapq', type=int, default=5, help='Minimum mapping quality')
    parser.add_argument('--min-len', type=int, default=50000, help='Minimum alignment length')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug output')
    
    args = parser.parse_args()
    
    liftover_bed_file(args.bed, args.paf, args.output, args.unmapped, args.debug)

if __name__ == '__main__':
    main()