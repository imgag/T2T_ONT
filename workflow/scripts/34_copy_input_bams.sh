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

    # Create the subfolder if it doesn't exist and copy the file to the new location only if it doesn't exist yet
    echo "if [ ! -f \"$subfolder/$new_filename\" ]; then mkdir -p \"$subfolder\" && cp \"$file\" \"$subfolder/$new_filename\"; fi"

done

# Path to the association table
ASSOC_TABLE="doc/tables/flowcell_biological_sample.tsv"

# Read the association table and group samples by their biological sample ID and processing system
declare -A samples
declare -A final_merged_samples
while IFS=$'\t' read -r name_ngsd name_external project_name run_flowcell_id processing_system_name; do
    # Skip header line and empty lines
    if [[ "$name_ngsd" == "name_ngsd" || -z "$name_ngsd" ]]; then
        continue
    fi

    if [[ $name_external == "T2T"* ]]; then
        # Create a composite key using biological sample and processing system
        composite_key="${name_external}_${processing_system_name}"
        if [[ -z "${samples[$composite_key]}" ]]; then
            samples["$composite_key"]="$name_ngsd"
        else
            samples["$composite_key"]="${samples[$composite_key]} $name_ngsd"
        fi
    fi
done < "$ASSOC_TABLE"

echo "pushd megSAP"

# Merge samples that belong to the same biological sample and have the same processing system
for composite_key in "${!samples[@]}"; do
    sample_list=(${samples["$composite_key"]})

    # Extract biological sample and processing system from the composite key
    biological_sample=$(echo "$composite_key" | sed 's/_UL-LR-ONT-SQK-ULK114$\|_PoreC-ONT-SQK-LSK114$\|_APK-ONT-SQK-APK114$\|_LR-ONT-SQK-LSK114$//')
    processing_system=$(echo "$composite_key" | sed -n 's/.*_\(UL-LR-ONT-SQK-ULK114\|PoreC-ONT-SQK-LSK114\|APK-ONT-SQK-APK114\|LR-ONT-SQK-LSK114\)$/\1/p')

    # Echo command for processing this biological sample with specific processing system
    echo "# Processing $biological_sample with processing system: $processing_system"

    # Sort the sample list to ensure sequential merging
    IFS=$'\n' sorted_samples=($(sort <<<"${sample_list[*]}"))
    unset IFS

    if [[ ${#sorted_samples[@]} -gt 1 ]]; then
        # Perform sequential merging
        for ((i=0; i<${#sorted_samples[@]}-1; i++)); do
            current_sample="${sorted_samples[$i]}"
            next_sample="${sorted_samples[$i+1]}"
            echo php /mnt/storage2/megSAP/pipeline/src/Tools/merge_samples.php -ps "$current_sample" -into "$next_sample"
        done

        # Store the last merged sample as the final one for this biological sample and processing system
        final_merged_sample="${sorted_samples[-1]}"
        if [[ $processing_system == "UL-LR-ONT-SQK-ULK114" ]]; then

            final_merged_samples["$biological_sample"]="$final_merged_sample"
        fi
    elif [[ ${#sorted_samples[@]} -eq 1 ]]; then
        # Only one sample for this biological sample and processing system
        if [[ $processing_system == "UL-LR-ONT-SQK-ULK114" ]]; then

            final_merged_samples["$biological_sample"]="${sorted_samples[0]}"
        fi
    fi
done



# Add command to run analysis on all merged UL-LR-ONT-SQK-ULK114 samples
echo
echo "# Run analysis on all merged UL-LR-ONT-SQK-ULK114 samples"


# If no samples were collected, print a warning
if [[ ${#final_merged_samples[@]} -eq 0 ]]; then
    echo "# WARNING: No UL-LR-ONT-SQK-ULK114 samples were collected"
else
    for biological_sample in "${!final_merged_samples[@]}"; do
        merged_sample="${final_merged_samples[$biological_sample]}"
        echo "php /mnt/storage2/users/ahgrosc1/dev/megsap/pipeline/src/Pipelines/analyze_longread.php -folder Sample_${merged_sample} -name ${merged_sample} -threads 32 # Processing biological sample: $biological_sample"
    done
fi
