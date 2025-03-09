# A Snakemake pipeline to run EvoTraceR and MACHINA analysis on data from CRISPR lineage tracing of cancer metastasis

This repo provides a snakemake pipeline to run the [EvoTraceR](https://github.com/Nowak-Lab/EvoTraceR) and [MACH2](https://github.com/elkebir-group/MACH2) analysis on input fastq files from the Nowak lab's data. To run, update necessary components in `config/config.yaml` and then run `./submit.sh` from the terminal.


### A note on the management of environments

I used a combination of singularity and conda to manage environments in the `Snakefile`. The singularity image files are too large to be uploaded, so they will have to be built or I can provide them on request. The step that runs MACH2 is not easy to setup in a way that works across platforms since MACH2 requires on a user specific Gurobi license. This step will have to be tuned to the user's setup.


This pipeline is not very easy to run across platforms in its current form and will be updated to make this easier.