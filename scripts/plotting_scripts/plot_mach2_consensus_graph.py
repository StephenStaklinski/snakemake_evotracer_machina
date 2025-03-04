#!/usr/bin/env python3

import sys
import networkx as nx


DEFAULT_COLORS = ["#006400", "#FF0000", "#0000CD", "#FFA500", "#800080", "#808080", "#FFC0CB", "#ADD8E6", "#A52A2A", "#FFFF00"]*3

# inputs
consensus_graph = sys.argv[1]
primary_tissue = sys.argv[2]
outfile = sys.argv[3]

# read the consensus graph
migrations = {}
all_tissues = []
with open(consensus_graph) as f:
    for line in f:
        migration, prob = line.strip().split(",")
        source, target, num = migration.split("_")
        all_tissues.append(source)
        all_tissues.append(target)
        migration_reduced = f"{source}_{target}"
        if migration_reduced in migrations:
            migrations[migration_reduced] += 1
        else:
            migrations[migration_reduced] = 1

# find all tissues to set the node colors
all_tissues = sorted(list(set(all_tissues) - {primary_tissue}))
custom_colors = {node: color for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)]) if node != primary_tissue}
all_tissues = [primary_tissue] + all_tissues
custom_colors[primary_tissue] = "black"

# plot the probability graph with edge thicknesses proportional to the probability
G = nx.MultiDiGraph()
for node in all_tissues:
    G.add_node(node, color=custom_colors[node], shape="box", fillcolor="white", penwidth=3.0, fontsize=32)
for edge, num in migrations.items():
    source, target = edge.split('_')
    label = ""
    if num > 1:
        label = f"{num}"
    G.add_edge(source, target, color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"', penwidth=3, label=label, fontsize=24)
dot = nx.nx_pydot.to_pydot(G)
dot.write_pdf(outfile)


