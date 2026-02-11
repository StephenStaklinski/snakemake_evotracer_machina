#!/bin/bash

dir="/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2"
outfile="/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2/all_mach2_data.csv"

consensusGraphs=$(find $dir -type f -name "consensus_graph.txt")

for file in $consensusGraphs; do
    mouse=$(basename $(dirname $(dirname $file)))
    cp=$(basename $(dirname $file))
    echo $file
    echo $mouse
    echo $cp

    while IFS= read -r line; do
        echo "$mouse,$cp,$line" | sed 's/_/,/g' >> "$outfile"
    done < "$file"
done
