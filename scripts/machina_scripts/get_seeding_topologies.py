#!/usr/bin/env python3

import sys
from ete3 import Tree

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
            metastatic_reseeding_children = set()
            
            # First check for metastatic reseeding events
            if tissue_dict[node.name] != ptissue:  
                # Track which reseeding events we've already counted
                counted_metastatic_reseeding_events = set()
                # Check if any child's tissue matches an ancestor's tissue
                for child in node.children:
                    if tissue_dict[child.name] != ptissue and tissue_dict[node.name] != tissue_dict[child.name]:
                        # Get all ancestor tissues
                        ancestor = node
                        while ancestor:
                            if tissue_dict[ancestor.name] == tissue_dict[child.name]:
                                # Create a unique identifier for this reseeding event
                                event_id = (tissue_dict[child.name], child.name)
                                if event_id not in counted_metastatic_reseeding_events:
                                    seeding_counts["metastatic_re_seeding"] += 1
                                    counted_metastatic_reseeding_events.add(event_id)
                                    metastatic_reseeding_children.add(child)
                                break
                            ancestor = ancestor.up

            # Categorize children, excluding those involved in metastatic reseeding
            for child in node.children:
                if child not in metastatic_reseeding_children:
                    if tissue_dict[child.name] != tissue_dict[node.name]:
                        children_diff.append(tissue_dict[child.name])
                    else:
                        children_same.append(tissue_dict[child.name])

            # If there are multiple different tissues among the children, it's a parallel seeding event
            if len(set(children_diff)) + len(metastatic_reseeding_children) > 1:
                if tissue_dict[node.name] == ptissue:
                    seeding_counts["primary_parallel_seeding"] += len(children_diff)
                # metastatic parallel seeding should not count primary reseeding event
                elif len(set(children_diff)) + len(metastatic_reseeding_children) > 2 or ptissue not in children_diff:
                    seeding_counts["metastatic_parallel_seeding"] += len([t for t in children_diff if t != ptissue])
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
        elif pair[0] != ptissue and pair[1] != ptissue and pair[0] != pair[1]:
            seeding_counts["metastatic_mono_seeding"] += 1
    
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
