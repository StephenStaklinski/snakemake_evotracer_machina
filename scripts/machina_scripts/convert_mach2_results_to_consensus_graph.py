#!/usr/bin/env python3

import sys

combined_graphs_file = sys.argv[1]
outfile = sys.argv[2]

migration_counts = {}
results_checked = set()
with open(combined_graphs_file, 'r') as f:
    for line in f.readlines():
        # skip the header line
        if "result_num" in line:
            continue
        result_num, source, target, count = line.strip().split(',')
        results_checked.add(result_num)
        for i in range(1, int(count) + 1):
            migration = f"{source}_{target}_{i}"
            if migration in migration_counts:
                migration_counts[migration] += 1
            else:
                migration_counts[migration] = 1
        
total_num_results = len(results_checked)

# normalize the counts to probabilities by dividing by the total number of results
normalized_migration_counts = {}
for migration in migration_counts:
    normalized_migration_counts[migration] = migration_counts[migration] / total_num_results

# write the consensus graph to a file
counts_outfile = outfile.replace('.csv', '_raw_counts.csv')
with open(counts_outfile, 'w') as f:
    for migration in migration_counts:
        f.write(f"{migration},{migration_counts[migration]}\n")

# write the normalized migration counts to a file
with open(outfile, 'w') as f:
    for migration in normalized_migration_counts:
        f.write(f"{migration},{normalized_migration_counts[migration]}\n")
