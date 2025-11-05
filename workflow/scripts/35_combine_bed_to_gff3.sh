#!/bin/bash

# Script to combine two BED files into a valid GFF3 file
# Usage: ./combine_bed_to_gff3.sh flagger.bed nucflag.bed output.gff3

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <flagger_bed> <nucflag_bed> <output_gff3>"
    echo ""
    echo "Arguments:"
    echo "  flagger_bed:  BED9 file with Flagger predictions (Err, Hap, Dup)"
    echo "  nucflag_bed:  BED4 file with Nucflag predictions (MISJOIN, HET)"
    echo "  output_gff3:  Output GFF3 file"
    exit 1
fi

FLAGGER_BED="$1"
NUCFLAG_BED="$2"
OUTPUT_GFF3="$3"

# Check input files exist
if [ ! -f "$FLAGGER_BED" ]; then
    echo "Error: Flagger BED file not found: $FLAGGER_BED"
    exit 1
fi

if [ ! -f "$NUCFLAG_BED" ]; then
    echo "Error: Nucflag BED file not found: $NUCFLAG_BED"
    exit 1
fi

echo "Combining BED files into GFF3 format..."
echo "  Flagger: $FLAGGER_BED"
echo "  Nucflag: $NUCFLAG_BED"
echo "  Output:  $OUTPUT_GFF3"

# Create GFF3 file with header
cat > "$OUTPUT_GFF3" << 'EOF'
##gff-version 3
##description: Combined assembly QC predictions from Flagger and Nucflag
##source: Flagger (Err, Hap, Dup) and Nucflag (MISJOIN, HET)
EOF

# Process Flagger BED file (BED9 format)
# Convert 0-based to 1-based, add source and type information
awk 'BEGIN {OFS="\t"}
/^track/ {next}
/^#/ {next}
NF >= 9 {
    chr = $1
    start = $2 + 1  # Convert to 1-based
    end = $3
    feature_type = $4
    score = ($5 != "" && $5 != "0") ? $5 : "."
    strand = ($6 != "" && $6 != ".") ? $6 : "."
    rgb = $9
    
    # Map feature types to SO terms (Sequence Ontology)
    if (feature_type == "Hap") {
        so_term = "haplotype_region"
        description = "Haploid region"
    } else if (feature_type == "Dup") {
        so_term = "duplication"
        description = "Duplicated region"
    } else if (feature_type == "Err") {
        so_term = "error_region"
        description = "Error region"
    } else {
        so_term = "region"
        description = feature_type
    }
    
    # Create unique ID
    id = sprintf("%s_flagger_%s_%d_%d", chr, feature_type, start, end)
    
    # Build attributes
    attributes = sprintf("ID=%s;Name=%s;source=Flagger;feature_type=%s;color=%s;description=%s", 
                        id, feature_type, feature_type, rgb, description)
    
    # Print GFF3 line
    print chr, "Flagger", so_term, start, end, score, strand, ".", attributes
}' "$FLAGGER_BED" >> "$OUTPUT_GFF3"

# Process Nucflag BED file (BED4 format)
awk 'BEGIN {OFS="\t"}
/^track/ {next}
/^#/ {next}
NF >= 4 {
    chr = $1
    start = $2 + 1  # Convert to 1-based
    end = $3
    feature_type = $4
    
    # Map feature types to SO terms
    if (feature_type == "MISJOIN") {
        so_term = "misassembly"
        description = "Potential misjoin"
        score = "."
    } else if (feature_type == "HET") {
        so_term = "heterozygous_region"
        description = "Heterozygous region"
        score = "."
    } else {
        so_term = "region"
        description = feature_type
        score = "."
    }
    
    # Create unique ID
    id = sprintf("%s_nucflag_%s_%d_%d", chr, feature_type, start, end)
    
    # Build attributes
    attributes = sprintf("ID=%s;Name=%s;source=Nucflag;feature_type=%s;description=%s", 
                        id, feature_type, feature_type, description)
    
    # Print GFF3 line
    print chr, "Nucflag", so_term, start, end, score, ".", ".", attributes
}' "$NUCFLAG_BED" >> "$OUTPUT_GFF3"

# Sort the GFF3 file by chromosome and position
echo "Sorting GFF3 file..."
(
    grep "^##" "$OUTPUT_GFF3"
    grep -v "^##" "$OUTPUT_GFF3" | sort -k1,1 -k4,4n
) > "${OUTPUT_GFF3}.tmp" && mv "${OUTPUT_GFF3}.tmp" "$OUTPUT_GFF3"

echo "Done! Created GFF3 file: $OUTPUT_GFF3"

# Print statistics
echo ""
echo "Statistics:"
echo "  Flagger features: $(grep -c "Flagger" "$OUTPUT_GFF3" || echo 0)"
echo "  Nucflag features: $(grep -c "Nucflag" "$OUTPUT_GFF3" || echo 0)"
echo "  Total features:   $(grep -cv "^##" "$OUTPUT_GFF3" || echo 0)"
