#!/usr/bin/env python3

import sys
import ete3
import pandas as pd

tree = sys.argv[1]
primary_tissue = str(sys.argv[2])
tissues_tsv = str(sys.argv[3])
output_dir = sys.argv[4]
outprefix = sys.argv[5]


# tree = "/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/billy_data_7_1_24_analysis_on_3_3_25/evotracer/MMUS1782/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick"
# primary_tissue = "BDR"
# tissues_tsv = "/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/billy_data_7_1_24_analysis_on_3_3_25/cp_split/MMUS1782/CP01/tissues.tsv"
# output_dir = "./"
# outprefix = "test"

# read in tree
tree = ete3.Tree(tree, format=0)

# read in tissues
tissue_data = pd.read_csv(tissues_tsv, sep="\t")

# prune tree to only include leaves with tissue data for this input CP
tips = tissue_data['group_name'].tolist()
tree.prune(tips)

# name internal nodes
i = 0
for node in tree.traverse():
    if node.is_root():
        node.name = "root"
    elif not node.is_leaf():
        node.name = f"node{i}"
        i += 1

# get edges and tip tissue labels expanded for multiple tissues
edges = []
tissue_labels = []
for node in tree.traverse():
    if node.is_root():
        # skip root since we look up for each edge and root has no parent to look up for
        continue
    elif not node.is_leaf():
        # add internal node edges
        edges.append((node.up.name, node.name))
    else:    

        leaf_tissues = tissue_data[tissue_data['group_name'] == node.name]['tissues'].values.tolist()[0]
        if "," in leaf_tissues:
            leaf_tissues = leaf_tissues.split(",")

            # remove multiple primaries
            leaf_tissues = list(set([primary_tissue if primary_tissue in tis else tis for tis in leaf_tissues]))

            # make sure there is more than 1 tissue left
            if len(leaf_tissues) == 1:
                edges.append((node.up.name, node.name))
                tissue_labels.append((node.name, leaf_tissues[0]))
            else:
                i = 0
                for tissue in leaf_tissues:
                    edges.append((node.up.name, f"{node.name}_{i}"))
                    tissue_labels.append((f"{node.name}_{i}", tissue))
                    i += 1
        else:
            # only one primary sample allowed here
            if primary_tissue in leaf_tissues:
                leaf_tissues = primary_tissue
            edges.append((node.up.name, node.name))
            tissue_labels.append((node.name, leaf_tissues))
            

tissues = list(set([x[1] for x in tissue_labels] + [primary_tissue]))
tissues.sort(key=lambda x: (x != primary_tissue, x))

i = 1
color_map = {}
for tissue in tissues:
    color_map[tissue] = i
    i += 1

# output files
output_file_leaf = output_dir + "/" + outprefix + ".labeling"
output_file_edges = output_dir + "/" + outprefix + ".tree"
output_file_colors = output_dir + "/" + outprefix + "_colors.txt"


with open(output_file_edges, "w") as file:
    for edge in edges:
        file.write(f'{edge[0]}\t{edge[1]}\n')

with open(output_file_leaf, "w") as file:
    for label in tissue_labels:
        file.write(f'{label[0]}\t{label[1]}\n')

with open(output_file_colors, "w") as file:
    for key, value in color_map.items():
        file.write(f'{key}\t{value}\n')


