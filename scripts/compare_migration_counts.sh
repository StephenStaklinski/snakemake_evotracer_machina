#!/bin/bash

indir="/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina/results/analysis_on_3_7_25_serio_pca_2p_2pr_data_3_7_25/mach2"

files=$(find $indir -type f -name "consensus_graph_filtered_by_threshold.txt")

outfile_counts="$indir/migration_counts.csv"
echo "mouse,cp,num_migrations" > $outfile_counts

for file in $files; do
    num_migrations=$(cat $file | wc -l)
    mouse=$(echo $file | awk -F'/' '{print $(NF-2)}')
    cp=$(echo $file | awk -F'/' '{print $(NF-1)}')
    echo "$mouse,$cp,$num_migrations" >> $outfile_counts
done

python scripts/plotting_scripts/plot_migration_counts.py $outfile_counts $indir/migration_counts.pdf

