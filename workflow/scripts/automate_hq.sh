for f in published.HQ*;do
	### mkdir -p "${f%.fastq.gz}"
	### mv $f "${f%.fastq.gz}"
	cd "${f%.fastq.gz}"
	rm bandage_info.txt
	conda run -n hifiasm hifiasm -t48 -o chr_19.asm --ont "$f".fastq.gz
	conda run -n bandage Bandage image chr_19.asm.bp.p_ctg.gfa main_contig_asm_graph.png
	conda run -n bandage Bandage info chr_19.asm.bp.p_ctg.gfa &> bandage_info.txt
	cd ..
done


for fol in published.HQ*;do
	hqname=$(echo $fol | sed 's/published.//g;s/.chr19//g')
	for f in published.UL*;do
		ulname=$(echo $f | sed 's/published.//g;s/.chr19.fastq.gz//g')
		folname=hybrid_"$hqname"_"$ulname"
		mkdir -p $folname
		cp $fol/*.fastq.gz $folname
		cp $f $folname
	done
done		

for f in hybrid*;do
	cd $f
	conda run -n hifiasm hifiasm -t48 -o chr_19.asm --ul published.UL* --ont published.HQ*
	conda run -n bandage Bandage image chr_19.asm.bp.p_ctg.gfa main_contig_asm_graph.png
	conda run -n bandage Bandage info chr_19.asm.bp.p_ctg.gfa &> bandage_info.txt
	rm *.bin
	cd ..
done

echo -e "Filename\tNode count\tEdge count\tSmallest edge overlap (bp)\tLargest edge overlap (bp)\tTotal length (bp)\tTotal length no overlaps (bp)\tDead ends\tPercentage dead ends\tConnected components\tLargest component (bp)\tTotal length orphaned nodes (bp)\tN50 (bp)\tShortest node (bp)\tLower quartile node (bp)\tMedian node (bp)\tUpper quartile node (bp)\tLongest node (bp)\tMedian depth\tEstimated sequence length (bp)" >> summary_"$proj".txt

for proj in published hybrid;do
	for f in "$proj"*; do
		if [ -d "$f" ]; then
			file=$f/bandage_info.txt
			if [ -f "$file" ]; then
				node_count=$(grep "Node count" "$file" | awk '{print $NF}')
				edge_count=$(grep "Edge count" "$file" | awk '{print $NF}')
				smallest_overlap=$(grep "Smallest edge overlap (bp)" "$file" | awk '{print $NF}')
				largest_overlap=$(grep "Largest edge overlap (bp)" "$file" | awk '{print $NF}')
				total_length=$(grep "Total length (bp)" "$file" | awk '{print $NF}')
				total_length_no_overlaps=$(grep "Total length no overlaps (bp)" "$file" | awk '{print $NF}')
				dead_ends=$(grep "Dead ends" "$file" | awk '{print $NF}')
				percentage_dead_ends=$(grep "Percentage dead ends" "$file" | awk '{print $NF}')
				connected_components=$(grep "Connected components" "$file" | awk '{print $NF}')
				largest_component=$(grep "Largest component (bp)" "$file" | awk '{print $NF}')
				orphaned_length=$(grep "Total length orphaned nodes (bp)" "$file" | awk '{print $NF}')
				n50=$(grep "N50 (bp)" "$file" | awk '{print $NF}')
				shortest_node=$(grep "Shortest node (bp)" "$file" | awk '{print $NF}')
				lower_quartile=$(grep "Lower quartile node (bp)" "$file" | awk '{print $NF}')
				median_node=$(grep "Median node (bp)" "$file" | awk '{print $NF}')
				upper_quartile=$(grep "Upper quartile node (bp)" "$file" | awk '{print $NF}')
				longest_node=$(grep "Longest node (bp)" "$file" | awk '{print $NF}')
				median_depth=$(grep "Median depth" "$file" | awk '{print $NF}')
				estimated_length=$(grep "Estimated sequence length (bp)" "$file" | awk '{print $NF}')

				# Print the extracted values for the current file, separated by tabs
				echo -e "$f\t$node_count\t$edge_count\t$smallest_overlap\t$largest_overlap\t$total_length\t$total_length_no_overlaps\t$dead_ends\t$percentage_dead_ends\t$connected_components\t$largest_component\t$orphaned_length\t$n50\t$shortest_node\t$lower_quartile\t$median_node\t$upper_quartile\t$longest_node\t$median_depth\t$estimated_length" >> summary_"$proj".txt
			else
				echo "Error: File not found: $file" >&2
			fi
		fi
	done
done
