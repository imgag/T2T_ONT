nucflag \
-i TUE_02_03UL.HQ_herro.50x.bam \
-d plot \
--overlay_regions ../repeatmasker/rm_summary/TUE_02_03UL_nucplot.bed \
--threads 12 \
--output_status TUE_02_03UL.nucflag_status.bed \
--output_misasm TUE_02_03UL.nucflag_misasm.bed \
> nucflag.log 2>&1


# Replace colons with underscores in all filenames in current directory
for file in plot/*; do
    if [[ -f "$file" ]]; then
        newname=$(echo "$file" | tr ':' '_')
        if [[ "$file" != "$newname" ]]; then
            mv "$file" "$newname"
        fi
    fi
done
``