import os, sys
from glob import glob

outDirs = []
for fastqDir in config['fastqDirs']:
    paths = glob('{}/*'.format(fastqDir))
    reads = [path for path in paths if os.path.isfile(path)]
    if len(reads) == 0:
        raise ValueError("No fastq files found in the input directory {fastqDir}.")
    outd = os.path.dirname(fastqDir)
    outDirs.append(outd)


rule all:
    input:
        expand("{outdir}/evotracer_output/evotracer.RData", outdir = outDirs),
        expand("{outdir}/evotracer_output/evotracer_graphs/hist_freq_indels.pdf", outdir = outDirs),
        expand("{outdir}/machina_output/all_results_extended.txt", outdir = outDirs),
        expand("{outdir}/machina_output/machina_graphs/machina_migration_plots/trans_mx_all.pdf", outdir = outDirs)

rule runEvotracer:
    input:
        evoFastqDir = "{outdir}/" + config['fastqDirPrefix']
    output:
        asvStat = "{outdir}/evotracer_output/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
        nwk = "{outdir}/evotracer_output/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
        rDataObject = "{outdir}/evotracer_output/evotracer.RData",
    params:
        evoOutDir = lambda wildcards: "{}/evotracer_output".format(wildcards.outdir)
    singularity:
        'envs/evotracer.sif'
    shell:
        """
        Rscript scripts/evotracer_scripts/evotracer.R {input.evoFastqDir} {params.evoOutDir}
        """

rule plotEvotracerResults:
    input:
        rDataObject = "{outdir}/evotracer_output/evotracer.RData",
    output:
        indels = "{outdir}/evotracer_output/evotracer_graphs/hist_freq_indels.pdf"
    params:
        evoPlotsOutDir = lambda wildcards: "{}/evotracer_output/evotracer_graphs".format(wildcards.outdir)
    singularity:
        "envs/evotracer_plotting.sif"
    shell:
        """
        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.1_hist_freq_indels.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.2_hist_freq_seq_length.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.3_hist_freq_site_marked.R {input.rDataObject} {params.evoPlotsOutDir};

        # Rscript scripts/plotting_scripts/3_clonal_population_analysis/05.1_cp_stat_vis.R {input.rDataObject} {params.evoPlotsOutDir};

        # Rscript scripts/plotting_scripts/4_phylogenetic_analysis/06.1_tree_msa_bubble_all_clones.R {input.rDataObject} {params.evoPlotsOutDir}
        """

rule runMachina:
    input:
        asvStat = "{outdir}/evotracer_output/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
        nwk = "{outdir}/evotracer_output/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
    output:
        extendedMachina = "{outdir}/machina_output/all_results_extended.txt",
        migrationMachina = "{outdir}/machina_output/migration.txt",
        seedingMachina = "{outdir}/machina_output/seeding_topology.txt"
    params:
        machinaScripts = "scripts/machina_scripts/",
        machinaOutPrefix = lambda wildcards: "{}/machina_output".format(wildcards.outdir),
        primaryTissue = config['primaryTissue'],
        batches = config['batches'],
        threads = config['threads']
    conda:
        "envs/machina.yaml"
    shell:
        """
        scripts/machina_scripts/run_machina.sh --infile {input.asvStat} --tree {input.nwk} --primary-tissue {params.primaryTissue} --scripts {params.machinaScripts} --prefix {params.machinaOutPrefix} --keep-first-cp --threads {params.threads} --batches {params.batches}
        """

rule plotMachinaResults:
    input:
        extendedMachina = "{outdir}/machina_output/all_results_extended.txt",
        migrationMachina = "{outdir}/machina_output/migration.txt",
        seedingMachina = "{outdir}/machina_output/seeding_topology.txt"
    output:
        machinaRateMatrix = "{outdir}/machina_output/machina_graphs/machina_migration_plots/trans_mx_all.pdf",
    params:
        machinaGraphsOutDir = lambda wildcards: "{}/machina_output/machina_graphs".format(wildcards.outdir)
    singularity:
        "envs/evotracer_plotting.sif"
    shell:
        """
        #Rscript scripts/plotting_scripts/5_machina_analysis/02_machina_tree_v1.R {input.extendedMachina} {params.machinaGraphsOutDir};

        Rscript scripts/plotting_scripts/5_machina_analysis/03_machina_migration_v1.R {input.migrationMachina} {params.machinaGraphsOutDir};

        Rscript scripts/plotting_scripts/5_machina_analysis/04_machina_seeding_topology_v1.R {input.seedingMachina} {params.machinaGraphsOutDir};
        """