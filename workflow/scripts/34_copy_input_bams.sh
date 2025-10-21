#!/bin/bash

# Define the source and destination directories
source_dir="data/basecalled/sup"
dest_dir="megSAP"

# Initialize a counter for the id
declare -A sample_counts

# Loop through each file in the source directory
for file in "$source_dir"/*/*.bam; do
    # Extract the sample name from the file name
    filename=$(basename "$file")
    sample=$(echo "$filename" | cut -d'_' -f1)

    # Increment the id for this sample
    sample_counts["$sample"]=$((sample_counts["$sample"] + 1))
    id=$(printf "%02d" ${sample_counts["$sample"]})

    # Define the subfolder and new file name
    subfolder="$dest_dir/Sample_${sample}_$id"
    new_filename="${sample}_$id.mod.unmapped.bam"

    # Create the subfolder if it doesn't exist and copy the file to the new location
    echo mkdir -p "$subfolder" "&&" cp "$file" "$subfolder/$new_filename"

done

# Path to the association table
ASSOC_TABLE="doc/tables/flowcell_biological_sample.tsv"

# Read the association table and group samples by their biological sample ID
declare -A samples
while IFS=$'\t' read -r name_ngsd name_external project_name run_flowcell_id; do
    if [[ $name_external == "GE-MED-T2T"* ]]; then
        if [[ -z "${samples[$name_external]}" ]]; then
            samples["$name_external"]="$name_ngsd"
        else
            samples["$name_external"]="${samples[$name_external]} $name_ngsd"
        fi
    fi
done < "$ASSOC_TABLE"

echo "pushd megSAP"

# Merge samples that belong to the same biological sample
for biological_sample in "${!samples[@]}"; do
    sample_list=(${samples["$biological_sample"]})

    # Sort the sample list to ensure sequential merging
    IFS=$'\n' sorted_samples=($(sort <<<"${sample_list[*]}"))
    unset IFS

    # Perform sequential merging
    for ((i=0; i<${#sorted_samples[@]}-1; i++)); do
        current_sample="${sorted_samples[$i]}"
        next_sample="${sorted_samples[$i+1]}"
        echo php /mnt/storage2/megSAP/pipeline/src/Tools/merge_samples.php -ps "$current_sample" -into "$next_sample"
    done
done
