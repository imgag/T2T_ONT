#!/usr/bin/env python3
"""
Liftover BED files using PAF alignments while preserving all original columns.
Enhanced version with better handling of unmappable regions.
Uses only the longest alignment for each region to avoid duplicates.
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

def find_nearest_mappable_position(pos, aln, direction='downstream', max_search=1000, debug=False):
    """
    Find the nearest mappable position when the exact position falls in insertion/deletion.
    """
    if not aln.cigar:
        # Simple linear mapping without CIGAR
        if aln.strand == '+':
            offset = pos - aln.query_start
            return aln.target_start + offset
        else:
            offset = aln.query_end - pos
            return aln.target_start + offset
    
    # Use CIGAR for precise mapping
    operations = parse_cigar_operations(aln.cigar)
    
    search_positions = []
    if direction == 'downstream':
        search_positions = [pos + i for i in range(max_search + 1)]
    else:  # upstream
        search_positions = [pos - i for i in range(max_search + 1) if pos - i >= aln.query_start]
    
    for search_pos in search_positions:
        if search_pos < aln.query_start or search_pos >= aln.query_end:
            continue
            
        query_pos = aln.query_start
        target_pos = aln.target_start
        
        for length, op in operations:
            if op == 'M':  # Match/mismatch
                if query_pos <= search_pos < query_pos + length:
                    if aln.strand == '+':
                        result = target_pos + (search_pos - query_pos)
                    else:
                        result = target_pos + length - (search_pos - query_pos) - 1
                    
                    if debug:
                        offset = abs(search_pos - pos)
                        print(f"    Found mappable position {offset}bp {direction} of original: {pos} -> {search_pos} -> {result}", file=sys.stderr)
                    return result
                query_pos += length
                target_pos += length
            elif op == 'I':  # Insertion in query
                query_pos += length
            elif op == 'D':  # Deletion in query
                target_pos += length
    
    if debug:
        print(f"    No mappable position found within {max_search}bp {direction} of {pos}", file=sys.stderr)
    return None

def liftover_position_robust(pos, aln, max_search=1000, debug=False):
    """
    Liftover a single position with fallback to nearest mappable position.
    """
    # First try exact mapping
    if not aln.cigar:
        # Simple linear mapping without CIGAR
        if aln.strand == '+':
            offset = pos - aln.query_start
            return aln.target_start + offset
        else:
            offset = aln.query_end - pos
            return aln.target_start + offset
    
    # Try exact CIGAR mapping first
    operations = parse_cigar_operations(aln.cigar)
    query_pos = aln.query_start
    target_pos = aln.target_start
    
    for length, op in operations:
        if op == 'M':  # Match/mismatch
            if query_pos <= pos < query_pos + length:
                if aln.strand == '+':
                    result = target_pos + (pos - query_pos)
                else:
                    result = target_pos + length - (pos - query_pos) - 1
                if debug:
                    print(f"    Exact mapping: {pos} -> {result}", file=sys.stderr)
                return result
            query_pos += length
            target_pos += length
        elif op == 'I':  # Insertion in query
            if query_pos <= pos < query_pos + length:
                # Position falls in insertion - try to find nearest mappable position
                if debug:
                    print(f"    Position {pos} falls in insertion, searching for nearest mappable position", file=sys.stderr)
                # Try downstream first (for gap regions, we want to map after the gap)
                result = find_nearest_mappable_position(pos, aln, 'downstream', max_search, debug)
                if result is not None:
                    return result
                # Try upstream as fallback
                return find_nearest_mappable_position(pos, aln, 'upstream', max_search, debug)
            query_pos += length
        elif op == 'D':  # Deletion in query
            target_pos += length
    
    # Position not found in any block - try nearest mappable position
    if debug:
        print(f"    Position {pos} not found in alignment, searching nearby", file=sys.stderr)
    result = find_nearest_mappable_position(pos, aln, 'downstream', max_search, debug)
    if result is not None:
        return result
    return find_nearest_mappable_position(pos, aln, 'upstream', max_search, debug)

def find_best_overlapping_alignment(bed_entry, alignments, window=10000):
    """
    Find the single best alignment that overlaps with the BED region.
    Prioritizes by: 1) Direct overlap, 2) Alignment length, 3) MAPQ
    
    Returns:
        (alignment, overlap_start, overlap_end, is_windowed) or None
    """
    if bed_entry.chrom not in alignments:
        return None
    
    best_alignment = None
    best_score = -1
    best_overlap = None
    best_windowed = False
    
    # First try direct overlaps
    for aln in alignments[bed_entry.chrom]:
        overlap_start = max(bed_entry.start, aln.query_start)
        overlap_end = min(bed_entry.end, aln.query_end)
        
        if overlap_start < overlap_end:
            # Calculate score: alignment_length * mapq * overlap_fraction
            overlap_length = overlap_end - overlap_start
            bed_length = bed_entry.end - bed_entry.start
            overlap_fraction = overlap_length / bed_length
            
            score = aln.alignment_len * aln.quality * overlap_fraction
            
            if score > best_score:
                best_score = score
                best_alignment = aln
                best_overlap = (overlap_start, overlap_end)
                best_windowed = False
    
    # If no direct overlaps found, try with window
    if best_alignment is None and window > 0:
        windowed_start = max(0, bed_entry.start - window)
        windowed_end = bed_entry.end + window
        
        for aln in alignments[bed_entry.chrom]:
            # Check if alignment overlaps with windowed region
            if aln.query_end > windowed_start and aln.query_start < windowed_end:
                # Calculate overlap with windowed region for scoring
                windowed_overlap_start = max(windowed_start, aln.query_start)
                windowed_overlap_end = min(windowed_end, aln.query_end)
                windowed_overlap_length = windowed_overlap_end - windowed_overlap_start
                
                # Score based on proximity and alignment quality
                distance_penalty = min(
                    abs(bed_entry.start - aln.query_end),
                    abs(bed_entry.end - aln.query_start)
                ) / window
                
                score = aln.alignment_len * aln.quality * (1 - distance_penalty)
                
                if score > best_score:
                    best_score = score
                    best_alignment = aln
                    
                    # Calculate actual mapping coordinates
                    actual_start = max(bed_entry.start, aln.query_start)
                    actual_end = min(bed_entry.end, aln.query_end)
                    
                    # If no direct overlap, use closest boundary
                    if actual_start >= actual_end:
                        if bed_entry.end <= aln.query_start:
                            # Region is upstream of alignment
                            actual_start = aln.query_start
                            actual_end = min(aln.query_start + (bed_entry.end - bed_entry.start), aln.query_end)
                        else:
                            # Region is downstream of alignment
                            actual_end = aln.query_end
                            actual_start = max(aln.query_end - (bed_entry.end - bed_entry.start), aln.query_start)
                    
                    best_overlap = (actual_start, actual_end)
                    best_windowed = True
    
    if best_alignment is not None:
        return (best_alignment, best_overlap[0], best_overlap[1], best_windowed)
    
    return None

def liftover_bed_region(bed_entry, alignments, max_search=1000, window=10000, debug=False):
    """
    Liftover a BED region using only the best overlapping alignment.
    """
    failure_reasons = []
    
    if bed_entry.chrom not in alignments:
        reason = f"No alignments found for chromosome '{bed_entry.chrom}'"
        failure_reasons.append(reason)
        if debug:
            print(f"  FAIL: {reason}", file=sys.stderr)
            print(f"        Available chromosomes: {list(alignments.keys())[:10]}...", file=sys.stderr)
        return [], failure_reasons
    
    # Find the single best overlapping alignment
    best_match = find_best_overlapping_alignment(bed_entry, alignments, window)
    
    if best_match is None:
        reason = f"No overlapping alignments found for region {bed_entry.start}-{bed_entry.end} even with {window}bp window"
        failure_reasons.append(reason)
        if debug:
            print(f"  FAIL: {reason}", file=sys.stderr)
        return [], failure_reasons
    
    aln, overlap_start, overlap_end, is_windowed = best_match
    
    if debug:
        windowed_str = " (windowed)" if is_windowed else ""
        print(f"    Best alignment: processing overlap {overlap_start}-{overlap_end}{windowed_str}", file=sys.stderr)
        print(f"      BED region: {bed_entry.start}-{bed_entry.end}", file=sys.stderr)
        print(f"      Alignment: {aln.query_start}-{aln.query_end} -> {aln.target_start}-{aln.target_end}", file=sys.stderr)
        print(f"      Alignment length: {aln.alignment_len}, MAPQ: {aln.quality}", file=sys.stderr)
    
    # Liftover start and end positions with robust mapping
    lifted_start = liftover_position_robust(overlap_start, aln, max_search, debug)
    lifted_end = liftover_position_robust(overlap_end - 1, aln, max_search, debug)  # BED end is exclusive
    
    if lifted_start is None:
        reason = f"Could not map start position {overlap_start} even with {max_search}bp search window"
        failure_reasons.append(reason)
        if debug:
            print(f"      FAIL: {reason}", file=sys.stderr)
        return [], failure_reasons
    
    if lifted_end is None:
        reason = f"Could not map end position {overlap_end-1} even with {max_search}bp search window"
        failure_reasons.append(reason)
        if debug:
            print(f"      FAIL: {reason}", file=sys.stderr)
        return [], failure_reasons
    
    # Ensure proper ordering
    if lifted_start > lifted_end:
        lifted_start, lifted_end = lifted_end, lifted_start
    
    lifted_end += 1  # Convert back to exclusive end
    
    if debug:
        print(f"      SUCCESS: {overlap_start}-{overlap_end} -> {lifted_start}-{lifted_end}", file=sys.stderr)
    
    # Create lifted BED entry preserving all original fields
    lifted_name = bed_entry.name
    if is_windowed:
        lifted_name += "_windowed"
    
    lifted_region = {
        'chrom': aln.target_name,
        'start': lifted_start,
        'end': lifted_end,
        'name': lifted_name,
        'score': bed_entry.score,
        'strand': aln.strand,
        'extra_fields': bed_entry.extra_fields,
        'original_entry': bed_entry,
        'is_windowed': is_windowed,
        'alignment_info': {
            'length': aln.alignment_len,
            'mapq': aln.quality,
            'target': f"{aln.target_name}:{aln.target_start}-{aln.target_end}"
        }
    }
    
    return [lifted_region], failure_reasons

def liftover_bed_file(bed_file, paf_file, output_file, unmapped_file, max_search=1000, window=10000, debug=False):
    """Main liftover function with enhanced parameters."""
    
    print(f"Loading PAF alignments from {paf_file}", file=sys.stderr)
    alignments = load_paf_alignments(paf_file)
    
    print(f"Found alignments for {len(alignments)} query sequences", file=sys.stderr)
    print(f"Using search window: {max_search}bp for unmappable positions", file=sys.stderr)
    print(f"Using alignment window: {window}bp for distant regions", file=sys.stderr)
    print(f"Strategy: Select single best alignment per region", file=sys.stderr)
    
    lifted_count = 0
    unmapped_count = 0
    windowed_count = 0
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
            
            if debug and line_num <= 10:  # Debug first 10 entries
                print(f"\nProcessing BED entry {line_num}: {bed_entry.chrom}:{bed_entry.start}-{bed_entry.end} ({bed_entry.name})", file=sys.stderr)
            
            lifted_regions, failure_reasons = liftover_bed_region(bed_entry, alignments, max_search, window, debug and line_num <= 10)
            
            if not lifted_regions:
                unmapped.write(line + '\n')
                
                # Write detailed failure reason to unmapped file
                if failure_reasons:
                    unmapped.write(f"# FAILURE REASONS: {'; '.join(failure_reasons)}\n")
                    for reason in failure_reasons:
                        failure_summary[reason] += 1
                
                unmapped_count += 1
                continue
            
            # Write lifted region (should be exactly one)
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
                
                if region['is_windowed']:
                    windowed_count += 1
    
    print(f"\nSUMMARY:", file=sys.stderr)
    print(f"Successfully lifted {lifted_count} regions", file=sys.stderr)
    print(f"  - Direct overlaps: {lifted_count - windowed_count}", file=sys.stderr)
    print(f"  - Using alignment window: {windowed_count}", file=sys.stderr)
    print(f"Failed to lift {unmapped_count} regions", file=sys.stderr)
    success_rate = (lifted_count / (lifted_count + unmapped_count)) * 100 if (lifted_count + unmapped_count) > 0 else 0
    print(f"Success rate: {success_rate:.1f}%", file=sys.stderr)
    
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
    parser.add_argument('--max-search', type=int, default=1000, help='Maximum distance to search for mappable positions (bp)')
    parser.add_argument('--window', type=int, default=10000, help='Window size for finding nearby alignments (bp)')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug output')
    
    args = parser.parse_args()
    
    liftover_bed_file(args.bed, args.paf, args.output, args.unmapped, args.max_search, args.window, args.debug)

if __name__ == '__main__':
    main()