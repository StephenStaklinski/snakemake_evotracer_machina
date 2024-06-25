# snakemake_evotracer_machina

This repo provides a snakemake pipeline to run the EvoTraceR to MACHINA analysis on input fastq files from Nowak Lab members. Please change the path to the fastqDirs in `config/config.yaml` and then run `./submit.sh` from the terminal.

### Managing environments

I used a combination of singularity and conda to manage environments in the `Snakefile`. The singularity image files are too large to be uploaded, so they will have to be built from the def files or I can provide them on request.

I made a singularity file for running evotracer with `envs/evotracer.def`. This retains the minimal installation to run evotracer.

I made a singularity file for plotting graphs with `envs/evotracer_plotting.def`.   

I uses conda with `envs/machina.yaml` in combination with CSHL HPC modules for running MACHINA, due to a dependency on an active gurobi license. Unfortunately, this will need to be specific to the users system, so if not using the CSHL HPC then one may need to follow the original [MACHINA instructions](https://github.com/raphael-group/machina) and integrate their own environment into the `rule runMachina:` in the `Snakefile`.