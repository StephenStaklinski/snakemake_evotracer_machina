# A Snakemake pipeline to run EvoTraceR and MACHINA analysis on data from CRISPR lineage tracing of cancer metastasis

This repo provides a snakemake pipeline to run the [EvoTraceR](https://github.com/Nowak-Lab/EvoTraceR) and [MACHINA](https://github.com/raphael-group/machina) analysis on input fastq files from the Nowak lab's BC10v0 style data. To run, change the path to the fastqDirs in `config/config.yaml` and then run `./submit.sh` from the terminal.

To run data from a different amplicon sequence, update `scripts/evotracer_scripts/evotracer.R`.  Further adjustments may be required in the `submit.sh` file depending on the cluster system used for submission, but this should be amenable to standard snakemake submission formats.

### A note on the management of environments

I used a combination of singularity and conda to manage environments in the `Snakefile`. The singularity image files are too large to be uploaded, so they will have to be built from the def files or I can provide them on request.

I made the singularity file `envs/evotracer.def` for running evotracer. This retains the minimal installation to run evotracer.

I made the singularity file `envs/evotracer_plotting.def` for plotting graphs.

I used the conda file `envs/machina.yaml` in combination with CSHL HPC modules for running MACHINA, due to a dependency on an active gurobi license. Unfortunately, this will need to be specific to the users system, so if not using the CSHL HPC then one may need to follow the original [MACHINA instructions](https://github.com/raphael-group/machina) and integrate their own environment into the `rule runMachina:` in the `Snakefile`.