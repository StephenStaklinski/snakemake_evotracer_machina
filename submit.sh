#!/bin/bash

export REPO_PATH=/grid/siepel/home_norepl/staklins/snakemake_evotracer_machina

snakemake \
-n \
--use-singularity \
--singularity-args "--bind $HOME/" \
--use-conda \
--snakefile $REPO_PATH/Snakefile \
--configfile $REPO_PATH/config/config.yaml \
--printshellcmds \
--keep-going \
--ignore-incomplete \
--cores 1 \
--jobs 1000 \
--cluster-config $REPO_PATH/config/cluster.yaml \
--cluster 'qsub -cwd -pe threads {cluster.cores} -l m_mem_free={cluster.mem} -o {cluster.logout} -e {cluster.logerror}'
