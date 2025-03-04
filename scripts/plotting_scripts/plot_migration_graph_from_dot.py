#!/usr/bin/env python3

import sys
import pydot


# Read command-line arguments
dot_file = sys.argv[1]
outfile = sys.argv[2]

# Read the DOT file
(graph,) = pydot.graph_from_dot_file(dot_file)

# Write the graph to a PDF file
graph.write_pdf(outfile)


