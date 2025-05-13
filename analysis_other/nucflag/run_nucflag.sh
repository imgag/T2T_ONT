awk -F'\t' 'BEGIN{OFS="\t"} $1!~/#/ {print $1,$4-1,$5,$3,$7}' analysis_other/repeatmasker/TUE_02_03UL/assembly.fasta.out.gff > repeatmasker_TUE_02_03UL.bed

nucflag -i TUE_02_03UL.HQ_herro.50x.bam -d plot --overlay_regions ../repeatmasker/TUE_02_03UL/assembly.fasta.out.bed
# Replace colons with underscores in all filenames in current directory
for file in plot/*; do
    if [[ -f "$file" ]]; then
        newname=$(echo "$file" | tr ':' '_')
        if [[ "$file" != "$newname" ]]; then
            mv "$file" "$newname"
        fi
    fi
done
