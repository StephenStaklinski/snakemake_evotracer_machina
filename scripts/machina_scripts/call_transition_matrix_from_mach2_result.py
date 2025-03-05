#!/usr/bin/env python3

import sys

treefile = sys.argv[1]
labelingfile = sys.argv[2]
outfile = sys.argv[3]


tissueMapping = {}
allTissues = set()
with open(labelingfile, "r") as f:
    for line in f:
        name, tissue = line.strip().split(' ')
        tissueMapping[name] = tissue
        allTissues.add(tissue)

transitionCounts = {tissue: {t: 0 for t in allTissues} for tissue in allTissues}
with open(treefile, "r") as f:
    for line in f:
        parent, child, num = line.strip().split(' ')
        parentTissue = tissueMapping[parent]
        childTissue = tissueMapping[child]
        transitionCounts[parentTissue][childTissue] += 1

# get the primary tissue from the input tree file name
primaryTissue = treefile.split("/")[-1].split("-")[0]

sortedTissues = sorted(allTissues - {primaryTissue})
sortedTissues.insert(0, primaryTissue)
sortedTransitionCounts = {parent: {child: transitionCounts[parent][child] for child in sortedTissues} for parent in sortedTissues}

header = "," + ",".join(list(sortedTransitionCounts.keys()))

with open(outfile, "w") as f:
    f.write(f"{header}\n")
    for parentTissue, children in sortedTransitionCounts.items():
        counts = f"{parentTissue}," + ",".join([str(count) for count in children.values()])
        f.write(f"{counts}\n")

