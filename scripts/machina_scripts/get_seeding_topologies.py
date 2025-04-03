#!/usr/bin/env python3

import sys
from collections import Counter
from ete3 import Tree
from itertools import groupby


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
    
    # Traverse the tree
    for i, node in enumerate(tree.traverse("preorder")):
        if not node.is_leaf():
            # Get the tissues of the children nodes and categorize them
            children_diff = []
            children_same = []
            for child in node.children:
                if tissue_dict[child.name] != tissue_dict[node.name]:
                    children_diff.append(tissue_dict[child.name])
                else:
                    children_same.append(tissue_dict[child.name])
            
            # If there are multiple different tissues among the children, it's a parallel seeding event
            if len(set(children_diff)) > 1:
                if tissue_dict[node.name] == ptissue:
                    seeding_counts["primary_parallel_seeding"] += 1
                # metastatic parallel seeding should not count primary reseeding event
                elif len(set(children_diff)) > 2 or ptissue not in children_diff:
                    seeding_counts["metastatic_parallel_seeding"] += 1
                # Add edges for only children that match the parent
                edges.extend([[tissue_dict[node.name], tissue_dict[node.name]] for _ in children_same])
            else:
                # Add edges for all children
                all_children = children_same + children_diff
                edges.extend([[tissue_dict[node.name], c] for c in all_children])

    # Process the edges to count different types of seeding events
    for pair in edges:
        if pair[0] == ptissue and pair[1] != ptissue:
            seeding_counts["primary_mono_seeding"] += 1
        elif pair[0] == ptissue and pair[1] == ptissue:
            seeding_counts["primary_confined"] += 1
        elif pair[0] != ptissue and pair[1] != ptissue and pair[0] == pair[1]:
            seeding_counts["metastatic_confined"] += 1
        elif pair[0] != ptissue and pair[1] == ptissue:
            seeding_counts["primary_re_seeding"] += 1
    
    # count metastatic mono seeding and metastatic reseeding events
    # ignore edges with primary tissue node and edges confined to one tissue
    filtered_edges = [sublist for sublist in edges if 'PRL' not in sublist]
    filtered_edges = [sublist for sublist in filtered_edges if len(set(sublist)) > 1]
    # get unique matching forward and reverse pairs
    forward_pairs = set()
    reverse_pairs = set()
    for pair in filtered_edges:
        tup = tuple(pair)
        rev = tuple(pair[::-1])
        if rev in forward_pairs:
            reverse_pairs.add(rev)
        forward_pairs.add(tup)
    
    # count metastatic reseeding events
    matching_pairs = [list(pair) for pair in reverse_pairs]
    for pair in matching_pairs:
        count = sum(1 for sublist in filtered_edges if sublist == pair)
        seeding_counts["metastatic_re_seeding"] += count
    
    # count remaining metastatic mono seeding events
    nonmatching_pairs = [sublist for sublist in filtered_edges if sublist not in matching_pairs]
    seeding_counts["metastatic_mono_seeding"] += len(nonmatching_pairs)

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
