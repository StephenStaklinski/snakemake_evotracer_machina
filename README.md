# Cancer metastasis CRISPR lineage tracing analysis pipeline

A Snakemake pipeline for analyzing CRISPR lineage tracing data from cancer metastasis studies using EvoTraceR and MACH2.

## Overview

This pipeline processes raw FASTQ files from CRISPR lineage tracing experiments to reconstruct and analyze tumor evolution and metastasis patterns. It combines several tools:

1. **EvoTraceR**: [EvoTraceR](https://github.com/Nowak-Lab/EvoTraceR) processes raw sequencing data and uses Cassiopeia-Greedy from [Cassiopeia](https://github.com/YosefLab/Cassiopeia) to reconstruct cell-lineage trees.
2. **MACH2**: [MACH2](https://github.com/elkebir-group/MACH2) infers a set of possible migration histories.
3. **BEAM**: [BEAM](https://github.com/CshlSiepelLab/BEAM) infers a joint cell-lineage tree and migration history posterior distribution using Markov Chain Monte Carlo (MCMC).
4. **VINE**: [VINE](https://github.com/CshlSiepelLab/VINE) approximates a joint cell-lineage tree and migration history posterior distribution using variational inference.

### Pipeline structure

The pipeline consists of several key components:

- `Snakefile`: Main pipeline definition
- `config/config.yaml`: Configuration file for pipeline parameters
- `envs/`: Environment definitions for different pipeline steps
- `scripts/`: Helper scripts for data processing


### Dependencies

The pipeline requires:
- Build singularity containers from within `envs/`
- Optain a Gurobi optimizer license (required for MACH2)
- Install BEAM

### Configuration

Before running the pipeline, update `config/config.yaml` with:
- Input FASTQ file paths
- Output directory
- Reference sequences and parameters
- Primary tissue definitions

## Usage

1. Configure the pipeline in `config/config.yaml`
2. Run the pipeline using:
   ```bash
   ./submit.sh
   ```

## Note on environment management

Due to size constraints, Singularity images are not included in the repository and must be built from the provided `.def` files or obtained separately (from staklins@cshl.edu). The MACH2 step requires a Gurobi license, which needs to be configured according to the user's setup. BEAM must be installed systemwide, as I have not yet made an effort to build BEAST2 in a container, which it requires.

Other changes may need to be made to the pipeline if not using it on the CSHL HPC server.