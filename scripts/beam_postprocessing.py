import sys
from graphposterior import PosteriorResults

trees_file = sys.argv[1]
state_key = sys.argv[2]
primary_tissue = sys.argv[3]
total_time = int(sys.argv[4])
burnin = float(sys.argv[5])
output_file_prefix = sys.argv[6]

# Initialize with BEAM output files
results = PosteriorResults(trees_file = trees_file, primary_tissue=primary_tissue, total_time=total_time, state_key=state_key, cores=1, burnin=burnin)

# Calculate consensus graph
results.get_consensus_graph(output_file=f"{output_file_prefix}_probability_graph.csv")

# Plot graphs
results.plot_thresholded_graph(threshold=[0.5, 0.75, 0.90], output_file_prefix=output_file_prefix)

