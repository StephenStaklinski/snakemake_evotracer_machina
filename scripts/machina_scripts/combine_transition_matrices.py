#!/use/bin/env python3

import sys

matrixFiles = sys.argv[1].split(',')
outfile = sys.argv[2]
normalized = bool(sys.argv[3])

overall_matrix = {}
all_tissues = set()
for file in matrixFiles:
    with open(file, 'r') as f:
        header_line = f.readline().strip()

        # skip empty files
        if header_line == '':
            continue

        recipient_tissue_order = header_line.split(',')[1:]
        for line in f:
            source_tissue, *counts = line.strip().split(',')
            if normalized:
                total_counts = sum(float(count) for count in counts)
                if total_counts == 0:
                    counts = [float(count) for count in counts]
                else:
                    counts = [float(count)/total_counts for count in counts]
            else:
                counts = [float(count) for count in counts]
            if source_tissue not in overall_matrix:
                overall_matrix[source_tissue] = {}
            for i, recipient_tissue in enumerate(recipient_tissue_order):
                if recipient_tissue not in overall_matrix[source_tissue]:
                    overall_matrix[source_tissue][recipient_tissue] = 0
                overall_matrix[source_tissue][recipient_tissue] += counts[i]
                all_tissues.add(recipient_tissue)
            all_tissues.add(source_tissue)

if normalized:
    for source_tissue in overall_matrix:
        recipient_sum = sum(overall_matrix[source_tissue].values())
        for recipient_tissue in overall_matrix[source_tissue]:
            overall_matrix[source_tissue][recipient_tissue] /= recipient_sum

# Ensure all tissues are in the matrix
sorted_tissues = sorted(list(all_tissues))
for tissue in sorted_tissues:
    if tissue not in overall_matrix:
        overall_matrix[tissue] = {}
    for recipient_tissue in sorted_tissues:
        if recipient_tissue not in overall_matrix[tissue]:
            overall_matrix[tissue][recipient_tissue] = 0

header = ',' + ','.join(sorted_tissues)

with open(outfile, 'w') as f:
    f.write(f'{header}' + '\n')
    for source_tissue in sorted_tissues:
        f.write(source_tissue + ',' + ','.join(str(overall_matrix[source_tissue][recipient_tissue]) for recipient_tissue in sorted_tissues) + '\n')

            