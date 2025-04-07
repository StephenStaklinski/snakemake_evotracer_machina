#!/usr/bin/env python3

import sys
from collections import Counter
from ete3 import Tree
from itertools import groupby


def find_cycles(edges, ptissue):
    """Find all cycles in a directed graph represented by edges, excluding any cycles that include the primary tissue."""
    # Create adjacency list
    graph = {}
    for src, dst in edges:
        if src not in graph:
            graph[src] = []
        graph[src].append(dst)
    
    # Track visited nodes and recursion stack
    visited = set()
    recursion_stack = set()
    cycles = []
    
    def dfs(node, path):
        visited.add(node)
        recursion_stack.add(node)
        path.append(node)
        
        if node in graph:
            for neighbor in graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, path)
                elif neighbor in recursion_stack:
                    # Found a cycle
                    cycle_start = path.index(neighbor)
                    cycle = path[cycle_start:]
                    # Only add cycle if it doesn't include primary tissue
                    if ptissue not in cycle:
                        cycles.append(cycle)
        
        path.pop()
        recursion_stack.remove(node)
    
    # Start DFS from each unvisited node
    for node in graph:
        if node not in visited:
            dfs(node, [])
    
    return cycles

def seeding_topology_tree(tabular_tree, tissue_dict, ptissue):
    # Convert the tabular tree into an ete3 Tree object
    tree = Tree.from_parent_child_table(tabular_tree)
    
    # Dictionary to count different types of seeding events
    seeding_counts = {
        "metastatic_confined": 0,
        "metastatic_mono_seeding": 0,
        "metastatic_parallel_seeding": 0,
        "metastatic_re_seeding": 0,
        "primary_confined": 0,
        "primary_mono_seeding": 0,
        "primary_parallel_seeding": 0,
        "primary_re_seeding": 0
    }

    # add the event from the origin to the root
    root = tree.get_tree_root().name
    if tissue_dict[root] == ptissue:
        seeding_counts["primary_confined"] += 1
    else:
        seeding_counts["primary_mono_seeding"] += 1
    
    # List to store edges for further processing
    edges = []
    
    for i, node in enumerate(tree.traverse("preorder")):
        if not node.is_leaf():
            for child in node.children:
                edges.append([tissue_dict[node.name], tissue_dict[child.name]])   
    
    # count metastatic mono seeding and metastatic reseeding events
    # ignore edges with primary tissue node and edges confined to one tissue
    filtered_edges = [sublist for sublist in edges if ptissue not in sublist]
    filtered_edges = [sublist for sublist in filtered_edges if len(set(sublist)) > 1]
    
    # Find all cycles in the filtered edges, explicitly excluding primary tissue
    cycles = find_cycles(filtered_edges, ptissue)
    
    # Count metastatic reseeding events based on cycles
    for cycle in cycles:
        # Each cycle represents a metastatic reseeding event
        seeding_counts["metastatic_re_seeding"] += 1
    
    # Create a set of edges that are part of cycles
    cycle_edges = set()
    for cycle in cycles:
        for i in range(len(cycle)):
            cycle_edges.add((cycle[i], cycle[(i+1)%len(cycle)]))
    
    # Count additional reseeding events between the same metastatic sites
    # First, count how many times each edge appears
    edge_counts = {}
    for edge in filtered_edges:
        edge_tuple = tuple(edge)
        edge_counts[edge_tuple] = edge_counts.get(edge_tuple, 0) + 1
    
    # For each edge that's part of a cycle, count additional occurrences as reseeding events
    for edge, count in edge_counts.items():
        if edge in cycle_edges and count > 1:
            seeding_counts["metastatic_re_seeding"] += (count - 1)
            
    # count metastatic parallel seeding events here
    non_cycle_edges = [edge for edge in filtered_edges if tuple(edge) not in cycle_edges]
    
    # Group edges by source tissue
    source_to_edges = {}
    for edge in non_cycle_edges:
        source = edge[0]
        if source not in source_to_edges:
            source_to_edges[source] = []
        source_to_edges[source].append(edge)
    
    # Count parallel seeding events
    for source, edges_by_src in source_to_edges.items():
        if len(edges_by_src) >= 2:
            # Get unique destination tissues
            dest_tissues = set(edge[1] for edge in edges_by_src)
            if len(dest_tissues) >= 2:
                seeding_counts["metastatic_parallel_seeding"] += len(edges_by_src)
            else:
                seeding_counts["metastatic_mono_seeding"] += len(edges_by_src)
        else:
            seeding_counts["metastatic_mono_seeding"] += len(edges_by_src)
    
    # count primary parallel and monoseeding events
    primary_edges = [edge for edge in edges if edge[0] == ptissue and edge[1] != ptissue]
    
    # Group primary edges by source tissue (should all be primary tissue)
    dest_tissues = set(edge[1] for edge in primary_edges)
    if len(dest_tissues) >= 2:
        seeding_counts["primary_parallel_seeding"] += len(primary_edges)
    else:
        seeding_counts["primary_mono_seeding"] += len(primary_edges)
    
    # Process the edges to count different types of seeding events
    for pair in edges:
        if pair[0] == ptissue and pair[1] == ptissue:
            seeding_counts["primary_confined"] += 1
        elif pair[0] != ptissue and pair[1] != ptissue and pair[0] == pair[1]:
            seeding_counts["metastatic_confined"] += 1
        elif pair[0] != ptissue and pair[1] == ptissue:
            seeding_counts["primary_re_seeding"] += 1

    # Return the counts of different seeding events
    return seeding_counts

infile_tree = sys.argv[1]
infile_tissues = sys.argv[2]
primary_tissue = sys.argv[3]
outfile = sys.argv[4]
cp = sys.argv[5]

# infile_tree = "/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/MMUS7317/CP10/PRL-T-0.tree"
# infile_tissues = "/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/MMUS7317/CP10/PRL-T-0.labeling"
# primary_tissue = "PRL"
# outfile = "/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/MMUS7317/CP10/PRL-T-0_seeding_topologies.txt"
# cp = "CP10"



tissueMapping = {}
allTissues = set()
allTissues.add(primary_tissue)
with open(infile_tissues, "r") as f:
    for line in f:
        try:
            name, tissue = line.strip().split('\t')
        except ValueError:
            name, tissue = line.strip().split(' ')
        tissueMapping[name] = tissue
        allTissues.add(tissue)

edges = []
with open(infile_tree, "r") as f:
    for line in f:
        try:
            parent, child, num = line.strip().split('\t')
        except ValueError:
            parent, child, num = line.strip().split(' ')
        edges.append((parent,child))

seedings = seeding_topology_tree(edges,tissueMapping,primary_tissue)

with open (outfile, "w") as f:
    f.write("CP,metastatic_confined,metastatic_mono_seeding,metastatic_parallel_seeding,metastatic_re_seeding,primary_confined,primary_mono_seeding,primary_parallel_seeding,primary_re_seeding\n")
    f.write(f"{cp},{seedings['metastatic_confined']},{seedings['metastatic_mono_seeding']},{seedings['metastatic_parallel_seeding']},{seedings['metastatic_re_seeding']},{seedings['primary_confined']},{seedings['primary_mono_seeding']},{seedings['primary_parallel_seeding']},{seedings['primary_re_seeding']}\n")
