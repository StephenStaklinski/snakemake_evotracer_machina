# Cancer metastasis CRISPR lineage tracing analysis pipeline

A Snakemake pipeline for analyzing CRISPR lineage tracing data from cancer metastasis studies using EvoTraceR and MACH2.

## Overview

This pipeline processes raw FASTQ files from CRISPR lineage tracing experiments to reconstruct and analyze tumor evolution and metastasis patterns. It combines two powerful tools:

1. **EvoTraceR**: [EvoTraceR](https://github.com/Nowak-Lab/EvoTraceR) processes raw sequencing data to reconstruct phylogenetic trees and analyze clonal evolution
2. **MACH2**: [MACH2](https://github.com/elkebir-group/MACH2) analyzes metastasis patterns and reconstructs seeding topologies

## Pipeline structure

The pipeline consists of several key components:

- `Snakefile`: Main pipeline definition
- `config/config.yaml`: Configuration file for pipeline parameters
- `envs/`: Environment definitions for different pipeline steps
- `scripts/`: Helper scripts for data processing

## Key features

- Automated processing of raw FASTQ files
- Phylogenetic tree reconstruction
- Clonal population analysis
- Metastasis pattern analysis
- Visualization of results including:
  - Clonal dispersal bargraphs
  - Consensus migration graphs
  - Transition matrices
  - Seeding topologies

## Setup and requirements

### Dependencies

The pipeline uses a combination of:
- Singularity containers for reproducible environments
- Conda environments for specific tools
- Gurobi optimizer (required for MACH2)

### Configuration

Before running the pipeline:
1. Update `config/config.yaml` with:
   - Input FASTQ file paths
   - Output directory
   - Reference sequences and parameters
   - Primary tissue definitions
2. Ensure you have the necessary Singularity images (contact staklins@cshl.edu if needed)
3. Set up Gurobi license for MACH2 analysis

## Usage

1. Configure the pipeline in `config/config.yaml`
2. Run the pipeline using:
   ```bash
   ./submit.sh
   ```

## Output

The pipeline generates several types of output files:
- Phylogenetic trees and character matrices
- Clonal population statistics
- Metastasis transition matrices
- Visualization plots in PDF format

## Note on environment management

The pipeline uses a combination of Singularity and Conda for environment management. Due to size constraints, Singularity images are not included in the repository and must be built or obtained separately. The MACH2 step requires a Gurobi license, which needs to be configured according to the user's setup.

Other changes may need to be made to the pipeline if not using it on the CSHL HPC server.