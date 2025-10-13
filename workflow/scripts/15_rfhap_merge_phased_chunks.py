#!/usr/bin/env python3

import argparse
import os
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


def parse_args():
    parser = argparse.ArgumentParser(description='Merge connected chunks with same haplotype assignment')
    parser.add_argument('-a', '--hap_a_list', required=True, help='File with paternal contig names')
    parser.add_argument('-b', '--hap_b_list', required=True, help='File with maternal contig names')
    parser.add_argument('-u', '--unassigned_list', required=True, help='File with unassigned contig names')
    parser.add_argument('-g', '--gfa', required=True, help='Original GFA file')
    parser.add_argument('-f', '--fasta', required=True, help='Original FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for files')
    parser.add_argument('--bridge_unassigned', action='store_true', 
                        help='Allow merging of haplotype chunks bridged by unassigned chunks')
    parser.add_argument('--max_unassigned_bridge', type=int, default=3,
                        help='Maximum number of consecutive unassigned chunks to use as a bridge (default: 3)')
    return parser.parse_args()


def load_haplotype_assignments(hap_a_file, hap_b_file, unassigned_file):
    """Load haplotype assignments from files"""
    haplotypes = {}
    
    with open(hap_a_file, 'r') as f:
        for line in f:
            chunk_id = line.strip()
            if chunk_id:
                haplotypes[chunk_id] = 'A'
    
    with open(hap_b_file, 'r') as f:
        for line in f:
            chunk_id = line.strip()
            if chunk_id:
                haplotypes[chunk_id] = 'B'
    
    with open(unassigned_file, 'r') as f:
        for line in f:
            chunk_id = line.strip()
            if chunk_id:
                haplotypes[chunk_id] = 'U'
    
    return haplotypes


def parse_gfa(gfa_file):
    """Parse GFA file into a graph and extract node information"""
    graph = nx.Graph()
    node_info = {}
    
    with open(gfa_file, 'r') as f:
        for line in f:
            if line.startswith('S'):  # Segment line
                fields = line.strip().split('\t')
                node_id = fields[1]
                
                # Extract contig origin and position info
                origin_contig = None
                start_pos = None
                end_pos = None
                
                for field in fields[3:]:
                    if field.startswith('OR:Z:'):
                        origin_contig = field[5:]
                    elif field.startswith('ST:i:'):
                        start_pos = int(field[5:])
                    elif field.startswith('EN:i:'):
                        end_pos = int(field[5:])
                
                node_info[node_id] = {
                    'origin_contig': origin_contig,
                    'start': start_pos,
                    'end': end_pos
                }
                graph.add_node(node_id)
                
            elif line.startswith('L'):  # Link line
                fields = line.strip().split('\t')
                from_node = fields[1]
                to_node = fields[3]
                graph.add_edge(from_node, to_node)
    
    return graph, node_info


def identify_merged_contigs(graph, node_info, haplotypes, bridge_unassigned=False, max_unassigned_bridge=3):
    """Identify connected components with the same haplotype, optionally bridging with unassigned chunks"""
    if not bridge_unassigned:
        # Standard approach: create subgraphs for each haplotype
        hap_a_nodes = [node for node, hap in haplotypes.items() if hap == 'A']
        hap_b_nodes = [node for node, hap in haplotypes.items() if hap == 'B']
        hap_u_nodes = [node for node, hap in haplotypes.items() if hap == 'U']
        
        hap_a_graph = graph.subgraph(hap_a_nodes)
        hap_b_graph = graph.subgraph(hap_b_nodes)
        hap_u_graph = graph.subgraph(hap_u_nodes)
        
        # Get connected components for each haplotype
        merged_contigs = {
            'A': list(nx.connected_components(hap_a_graph)),
            'B': list(nx.connected_components(hap_b_graph)),
            'U': list(nx.connected_components(hap_u_graph))
        }
    else:
        # Advanced approach: allow bridging through unassigned chunks
        merged_contigs = {'A': [], 'B': [], 'U': []}
        
        # Create a working copy of the graph
        working_graph = graph.copy()
        
        # Process each haplotype (A and B)
        for hap in ['A', 'B']:
            # Get nodes for this haplotype
            hap_nodes = [node for node, h in haplotypes.items() if h == hap]
            
            # For each node in this haplotype
            processed_nodes = set()
            for start_node in hap_nodes:
                if start_node in processed_nodes:
                    continue
                
                # Start a new component
                component = {start_node}
                processed_nodes.add(start_node)
                
                # Find all connected nodes of the same haplotype, possibly through unassigned bridges
                queue = [start_node]
                while queue:
                    current = queue.pop(0)
                    
                    # Check all neighbors up to max_unassigned_bridge steps away
                    for path_length in range(1, max_unassigned_bridge + 2):
                        # Find all nodes exactly path_length steps away
                        if path_length == 1:
                            neighbors = list(working_graph.neighbors(current))
                        else:
                            # Use networkx to find nodes at exact distance
                            neighbors = [n for n in nx.single_source_shortest_path_length(working_graph, current, cutoff=path_length) 
                                        if nx.single_source_shortest_path_length(working_graph, current, cutoff=path_length)[n] == path_length]
                        
                        for neighbor in neighbors:
                            # If it's the same haplotype and not processed
                            if haplotypes.get(neighbor) == hap and neighbor not in processed_nodes:
                                component.add(neighbor)
                                processed_nodes.add(neighbor)
                                queue.append(neighbor)
                                
                                # Check if there's a path through unassigned chunks
                                if path_length > 1:
                                    # Get the path
                                    path = nx.shortest_path(working_graph, current, neighbor)
                                    
                                    # Check if intermediate nodes are unassigned
                                    all_unassigned = all(haplotypes.get(node) == 'U' for node in path[1:-1])
                                    
                                    if all_unassigned:
                                        # Add the unassigned nodes to this component
                                        for bridge_node in path[1:-1]:
                                            component.add(bridge_node)
                                            # Mark as processed so it's not used in the 'U' haplotype
                                            processed_nodes.add(bridge_node)
                
                # Add the component to the results
                if component:
                    merged_contigs[hap].append(component)
            
            # Process remaining unassigned nodes
            unassigned_nodes = [node for node, h in haplotypes.items() if h == 'U' and node not in processed_nodes]
            unassigned_graph = working_graph.subgraph(unassigned_nodes)
            
            for component in nx.connected_components(unassigned_graph):
                if component:
                    merged_contigs['U'].append(component)
    
    return merged_contigs


def group_chunks_by_original_contig(node_info):
    """Group chunks by their original contig"""
    original_contigs = {}
    
    for chunk_id, info in node_info.items():
        origin = info['origin_contig']
        if origin not in original_contigs:
            original_contigs[origin] = []
        original_contigs[origin].append(chunk_id)
    
    # Sort chunks within each original contig by position
    for origin, chunks in original_contigs.items():
        original_contigs[origin] = sorted(chunks, key=lambda x: node_info[x]['start'])
    
    return original_contigs


def sort_chunks_by_position(chunks, node_info):
    """Sort chunks by their position in the original contig"""
    sorted_chunks = []
    
    # Group chunks by original contig
    contigs = {}
    for chunk in chunks:
        origin = node_info[chunk]['origin_contig']
        if origin not in contigs:
            contigs[origin] = []
        contigs[origin].append(chunk)
    
    # Sort chunks within each contig by start position
    for contig, contig_chunks in contigs.items():
        sorted_contig_chunks = sorted(contig_chunks, key=lambda x: node_info[x]['start'])
        sorted_chunks.extend(sorted_contig_chunks)
    
    return sorted_chunks


def create_merged_fasta(merged_contigs, node_info, fasta_file, output_file):
    """Create a FASTA file with merged contigs"""
    # Load original sequences
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    
    # Create merged sequences
    merged_records = []
    
    for hap in ['A', 'B', 'U']:
        for i, component in enumerate(merged_contigs[hap]):
            component_list = list(component)
            sorted_chunks = sort_chunks_by_position(component_list, node_info)
            
            # Get original contigs and positions
            contig_info = []
            for chunk in sorted_chunks:
                origin = node_info[chunk]['origin_contig']
                start = node_info[chunk]['start']
                end = node_info[chunk]['end']
                contig_info.append(f"{origin}:{start}-{end}")
            
            # Create merged sequence
            merged_seq = ""
            for chunk in sorted_chunks:
                merged_seq += sequences[chunk]
            
            # Create record
            merged_id = f"merged_hap{hap}_{i+1}"
            description = f"Merged from chunks: {','.join(sorted_chunks)} | Original contigs: {','.join(contig_info)}"
            record = SeqRecord(
                Seq(merged_seq),
                id=merged_id,
                description=description
            )
            merged_records.append(record)
    
    # Write to file
    SeqIO.write(merged_records, output_file, "fasta")


def get_bandage_color_tag(hap):
    """Get color tag for Bandage visualization"""
    if hap == 'A':
        return "CL:Z:blue"  # Paternal
    elif hap == 'B':
        return "CL:Z:red"   # Maternal
    else:
        return "CL:Z:gray"  # Unassigned


def get_haplotype_name(hap):
    """Get human-readable haplotype name"""
    if hap == 'A':
        return "Pat"  # Paternal
    elif hap == 'B':
        return "Mat"  # Maternal
    else:
        return "Unk"  # Unassigned/Unknown


def create_merged_gfa(merged_contigs, node_info, output_file):
    """Create a GFA file with merged contigs"""
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments
        for hap in ['A', 'B', 'U']:
            for i, component in enumerate(merged_contigs[hap]):
                component_list = list(component)
                sorted_chunks = sort_chunks_by_position(component_list, node_info)
                
                # Calculate length of merged segment
                length = sum(node_info[chunk]['end'] - node_info[chunk]['start'] + 1 for chunk in sorted_chunks)
                
                # Create segment ID and write
                segment_id = f"merged_hap{hap}_{i+1}"
                color_tag = get_bandage_color_tag(hap)
                hap_name = get_haplotype_name(hap)
                
                f.write(f"S\t{segment_id}\t*\tLN:i:{length}\tHP:Z:{hap}\t" +
                        f"HN:Z:{hap_name}\tOC:Z:{','.join(sorted_chunks)}\t{color_tag}\n")
        
        # No links in the merged GFA as each component becomes a single contig


def create_merged_noseq_gfa(merged_contigs, node_info, output_file):
    """Create a GFA file without sequences for merged contigs"""
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments (same as regular GFA but explicitly without sequences)
        for hap in ['A', 'B', 'U']:
            for i, component in enumerate(merged_contigs[hap]):
                component_list = list(component)
                sorted_chunks = sort_chunks_by_position(component_list, node_info)
                
                # Calculate length of merged segment
                length = sum(node_info[chunk]['end'] - node_info[chunk]['start'] + 1 for chunk in sorted_chunks)
                
                # Create segment ID and write
                segment_id = f"merged_hap{hap}_{i+1}"
                color_tag = get_bandage_color_tag(hap)
                hap_name = get_haplotype_name(hap)
                
                f.write(f"S\t{segment_id}\t*\tLN:i:{length}\tHP:Z:{hap}\t" +
                        f"HN:Z:{hap_name}\tOC:Z:{','.join(sorted_chunks)}\t{color_tag}\n")


def create_comprehensive_gfa(merged_contigs, node_info, haplotypes, output_file):
    """Create a comprehensive GFA file with all connections between chunks from the same original contig"""
    # Create a mapping from chunk_id to merged_id
    chunk_to_merged = {}
    for hap in ['A', 'B', 'U']:
        for i, component in enumerate(merged_contigs[hap]):
            merged_id = f"merged_hap{hap}_{i+1}"
            for chunk in component:
                chunk_to_merged[chunk] = merged_id
    
    # Group chunks by original contig
    original_contigs = group_chunks_by_original_contig(node_info)
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments (merged contigs)
        merged_info = {}
        for hap in ['A', 'B', 'U']:
            for i, component in enumerate(merged_contigs[hap]):
                component_list = list(component)
                sorted_chunks = sort_chunks_by_position(component_list, node_info)
                
                # Calculate length of merged segment
                length = sum(node_info[chunk]['end'] - node_info[chunk]['start'] + 1 for chunk in sorted_chunks)
                
                # Create segment ID and write
                segment_id = f"merged_hap{hap}_{i+1}"
                merged_info[segment_id] = {
                    'chunks': sorted_chunks,
                    'haplotype': hap,
                    'length': length
                }
                
                # Get original contigs
                original_contig_ids = set()
                for chunk in sorted_chunks:
                    original_contig_ids.add(node_info[chunk]['origin_contig'])
                
                # Add color and haplotype name tags for Bandage
                color_tag = get_bandage_color_tag(hap)
                hap_name = get_haplotype_name(hap)
                
                # Write segment line
                f.write(f"S\t{segment_id}\t*\tLN:i:{length}\tHP:Z:{hap}\t" +
                        f"HN:Z:{hap_name}\tOC:Z:{','.join(sorted_chunks)}\t" +
                        f"OR:Z:{','.join(original_contig_ids)}\t{color_tag}\n")
        
        # Write links between merged segments from the same original contig
        link_id = 1
        for origin, chunks in original_contigs.items():
            # Sort chunks by position
            sorted_chunks = sorted(chunks, key=lambda x: node_info[x]['start'])
            
            # Find consecutive chunks with different merged IDs
            for i in range(len(sorted_chunks) - 1):
                curr_chunk = sorted_chunks[i]
                next_chunk = sorted_chunks[i + 1]
                
                if curr_chunk in chunk_to_merged and next_chunk in chunk_to_merged:
                    curr_merged = chunk_to_merged[curr_chunk]
                    next_merged = chunk_to_merged[next_chunk]
                    
                    # Only add link if they are different merged segments
                    if curr_merged != next_merged:
                        curr_hap = merged_info[curr_merged]['haplotype']
                        next_hap = merged_info[next_merged]['haplotype']
                        
                        # Write link line
                        f.write(f"L\t{curr_merged}\t+\t{next_merged}\t+\t0M\t" +
                                f"LT:Z:haplotype_switch\tLH:Z:{curr_hap}->{next_hap}\tOC:Z:{origin}\tID:i:{link_id}\n")
                        link_id += 1


def create_comprehensive_noseq_gfa(merged_contigs, node_info, haplotypes, output_file):
    """Create a comprehensive GFA file without sequences"""
    # Just reuse the comprehensive GFA function since we're not including sequences anyway
    create_comprehensive_gfa(merged_contigs, node_info, haplotypes, output_file)


def create_merge_report(merged_contigs, node_info, output_file):
    """Create a report with information about merged contigs"""
    with open(output_file, 'w') as f:
        f.write("merged_id\thaplotype\thaplotype_name\tnum_chunks\toriginal_contigs\ttotal_length\n")
        
        for hap in ['A', 'B', 'U']:
            for i, component in enumerate(merged_contigs[hap]):
                component_list = list(component)
                sorted_chunks = sort_chunks_by_position(component_list, node_info)
                
                # Get original contigs
                original_contigs = set()
                for chunk in sorted_chunks:
                    original_contigs.add(node_info[chunk]['origin_contig'])
                
                # Calculate total length
                total_length = sum(node_info[chunk]['end'] - node_info[chunk]['start'] + 1 for chunk in sorted_chunks)
                
                # Get human-readable haplotype name
                hap_name = get_haplotype_name(hap)
                
                # Write report line
                merged_id = f"merged_hap{hap}_{i+1}"
                f.write(f"{merged_id}\t{hap}\t{hap_name}\t{len(sorted_chunks)}\t{','.join(original_contigs)}\t{total_length}\n")


def create_haplotype_switch_report(merged_contigs, node_info, haplotypes, output_file):
    """Create a report showing haplotype switches within original contigs"""
    # Group chunks by original contig
    original_contigs = group_chunks_by_original_contig(node_info)
    
    # Create a mapping from chunk_id to merged_id
    chunk_to_merged = {}
    for hap in ['A', 'B', 'U']:
        for i, component in enumerate(merged_contigs[hap]):
            merged_id = f"merged_hap{hap}_{i+1}"
            for chunk in component:
                chunk_to_merged[chunk] = merged_id
    
    with open(output_file, 'w') as f:
        f.write("original_contig\tnum_chunks\thaplotype_switches\tswitch_positions\thaplotype_sequence\thaplotype_names\n")
        
        for origin, chunks in original_contigs.items():
            # Sort chunks by position
            sorted_chunks = sorted(chunks, key=lambda x: node_info[x]['start'])
            
            # Track haplotype switches
            haplotype_sequence = []
            haplotype_names = []
            switch_positions = []
            
            for i, chunk in enumerate(sorted_chunks):
                hap = haplotypes.get(chunk, 'Unknown')
                haplotype_sequence.append(hap)
                haplotype_names.append(get_haplotype_name(hap))
                
                # Check for switch
                if i > 0 and haplotype_sequence[i] != haplotype_sequence[i-1]:
                    prev_pos = node_info[sorted_chunks[i-1]]['end']
                    curr_pos = node_info[chunk]['start']
                    switch_positions.append(f"{prev_pos}-{curr_pos}")
            
            # Count switches
            switches = 0
            for i in range(1, len(haplotype_sequence)):
                if haplotype_sequence[i] != haplotype_sequence[i-1]:
                    switches += 1
            
            # Write report line
            f.write(f"{origin}\t{len(sorted_chunks)}\t{switches}\t" +
                    f"{','.join(switch_positions)}\t{''.join(haplotype_sequence)}\t" +
                    f"{','.join(haplotype_names)}\n")


def create_bridging_report(merged_contigs, haplotypes, output_file):
    """Create a report showing which unassigned chunks were used as bridges"""
    with open(output_file, 'w') as f:
        f.write("merged_id\thaplotype\tnum_chunks\tnum_bridges\tbridge_chunks\n")
        
        for hap in ['A', 'B']:
            for i, component in enumerate(merged_contigs[hap]):
                component_list = list(component)
                
                # Find unassigned chunks used as bridges
                bridge_chunks = [chunk for chunk in component_list if haplotypes.get(chunk) == 'U']
                
                # Write report line
                merged_id = f"merged_hap{hap}_{i+1}"
                f.write(f"{merged_id}\t{hap}\t{len(component_list)}\t{len(bridge_chunks)}\t{','.join(bridge_chunks)}\n")


def main():
    args = parse_args()
    
    # Load haplotype assignments
    haplotypes = load_haplotype_assignments(args.hap_a_list, args.hap_b_list, args.unassigned_list)
    
    # Parse GFA file
    graph, node_info = parse_gfa(args.gfa)
    
    # Identify merged contigs
    merged_contigs = identify_merged_contigs(graph, node_info, haplotypes, 
                                            bridge_unassigned=args.bridge_unassigned,
                                            max_unassigned_bridge=args.max_unassigned_bridge)
    
    # Create merged FASTA
    merged_fasta = f"{args.output}.merged.fasta"
    create_merged_fasta(merged_contigs, node_info, args.fasta, merged_fasta)
    
    # Create merged GFA
    merged_gfa = f"{args.output}.merged.gfa"
    create_merged_gfa(merged_contigs, node_info, merged_gfa)
    
    # Create merged GFA without sequences
    merged_noseq_gfa = f"{args.output}.merged.noseq.gfa"
    create_merged_noseq_gfa(merged_contigs, node_info, merged_noseq_gfa)
    
    # Create comprehensive GFA with connections between original contigs
    comprehensive_gfa = f"{args.output}.comprehensive.gfa"
    create_comprehensive_gfa(merged_contigs, node_info, haplotypes, comprehensive_gfa)
    
    # Create comprehensive GFA without sequences
    comprehensive_noseq_gfa = f"{args.output}.comprehensive.noseq.gfa"
    create_comprehensive_noseq_gfa(merged_contigs, node_info, haplotypes, comprehensive_noseq_gfa)
    
    # Create merge report
    merge_report = f"{args.output}.merge_report.tsv"
    create_merge_report(merged_contigs, node_info, merge_report)
    
    # Create haplotype switch report
    switch_report = f"{args.output}.switch_report.tsv"
    create_haplotype_switch_report(merged_contigs, node_info, haplotypes, switch_report)
    
    # Create bridging report if bridging was enabled
    if args.bridge_unassigned:
        bridge_report = f"{args.output}.bridge_report.tsv"
        create_bridging_report(merged_contigs, haplotypes, bridge_report)
        print(f"Bridge report written to: {bridge_report}")
    
    # Print summary
    total_merged = sum(len(comps) for comps in merged_contigs.values())
    print(f"Merged chunks into {total_merged} contigs:")
    print(f"  - Haplotype A (Pat): {len(merged_contigs['A'])} contigs")
    print(f"  - Haplotype B (Mat): {len(merged_contigs['B'])} contigs")
    print(f"  - Unassigned: {len(merged_contigs['U'])} contigs")
    
    if args.bridge_unassigned:
        print(f"Bridging through unassigned chunks: ENABLED (max {args.max_unassigned_bridge} chunks)")
    else:
        print("Bridging through unassigned chunks: DISABLED")
        
    print(f"Merged FASTA written to: {merged_fasta}")
    print(f"Merged GFA written to: {merged_gfa}")
    print(f"Merged GFA without sequences written to: {merged_noseq_gfa}")
    print(f"Comprehensive GFA with connections written to: {comprehensive_gfa}")
    print(f"Comprehensive GFA without sequences written to: {comprehensive_noseq_gfa}")
    print(f"Merge report written to: {merge_report}")
    print(f"Haplotype switch report written to: {switch_report}")


if __name__ == "__main__":
    main() 