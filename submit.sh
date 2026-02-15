#!/bin/bash

snakemake \
-n \
--use-singularity \
--singularity-args "--bind $HOME/ --bind $GRB_LICENSE_FILE:/mnt/gurobi.lic --env GRB_LICENSE_FILE=/mnt/gurobi.lic" \
--snakefile ./Snakefile \
--configfile ./config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 1 \
--jobs 10000 \
--cluster-config ./config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -o {cluster.logout} -e {cluster.logerror}'
