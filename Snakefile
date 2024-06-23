import os
from glob import glob
from numpy import unique

paths = glob('{}/*'.format(config['fastqDir']))
reads = [path for path in paths if os.path.isfile(path)]
if len(reads) == 0:
    raise ValueError("No fastq files found in the input directory.")
mouseID = os.path.basename(reads[0]).split("_")[0]
tissues = unique([os.path.basename(path).split("_")[1] for path in reads])
sample = os.path.basename(reads[0]).split("_")[-2]
outdir = os.path.dirname(config['fastqDir'])

# # output extracted info from input fastq dir
# print("Info extracted from input fastq dir:")
# print(f"Fastq files: {reads}")
# print(f"Mouse: {mouseID}")
# print(f"Tissues: {tissues}")
# print(f"Sample: {sample}")
# print(f"Outdir: {outdir}")

rule all:
    input:
        f"{outdir}/evotracer_output/{mouseID}_evotracer.RData",
        f"{outdir}/evotracer_output/evotracer_graphs/cp_tree_msa_cna_bc_bubble_qnt_ggtree_mp.pdf"
        # f"{outdir}/machina_output/{mouseID}_all_results_extended.txt",
        # "test.txt"

rule runEvotracer:
    input:
        expand(config['fastqDir'] + '/{id}_{tissue}_BC10v0_MG_{sample}_R1.fastq', id=mouseID, tissue=tissues, sample=sample),
        expand(config['fastqDir'] + '/{id}_{tissue}_BC10v0_MG_{sample}_R2.fastq', id=mouseID, tissue=tissues, sample=sample)
    output:
        asvStat = outdir + "/evotracer_output/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
        nwk = outdir + "/evotracer_output/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
        rDataObject = f"{outdir}/evotracer_output/{mouseID}_evotracer.RData"
    params:
        fastqDir = config['fastqDir'],
        evoOutDir = outdir + "/evotracer_output"
    singularity:
        'envs/evotracer.sif'
    shell:
        """
        Rscript scripts/evotracer_scripts/evotracer.R {params.fastqDir} {params.evoOutDir} {mouseID}
        """

rule plotEvotracerResults:
    input:
        rDataObject = f"{outdir}/evotracer_output/{mouseID}_evotracer.RData"
    output:
        indels = f"{outdir}/evotracer_output/evotracer_graphs/hist_freq_indels.pdf",
        seqLen = f"{outdir}/evotracer_output/evotracer_graphs/hist_freq_seq_length.pdf",
        histFreq = f"{outdir}/evotracer_output/evotracer_graphs/hist_freq_site_affected.pdf",
        cpStats = f"{outdir}/evotracer_output/evotracer_graphs/stat_cps_dispersal_bargraph_hm.pdf",
        cpTree = f"{outdir}/evotracer_output/evotracer_graphs/cp_tree_msa_cna_bc_bubble_qnt_ggtree_mp.pdf"
    params:
        evoOutDir = outdir + "/evotracer_output/"
    singularity:
        "envs/evotracer_plotting.sif"
    shell:
        """
        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.1_hist_freq_indels.R {input.rDataObject} {params.evoOutDir};

        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.2_hist_freq_seq_length.R {input.rDataObject} {params.evoOutDir};

        Rscript scripts/plotting_scripts/2_barcode_edits_analysis/04.3_hist_freq_site_marked.R {input.rDataObject} {params.evoOutDir};

        Rscript scripts/plotting_scripts/3_clonal_population_analysis/05.1_cp_stat_vis.R {input.rDataObject} {params.evoOutDir};

        Rscript scripts/plotting_scripts/4_phylogenetic_analysis/06.1_tree_msa_bubble_all_clones.R {input.rDataObject} {params.evoOutDir}
        """

# rule runMachina:
#     input:
#         asvStat = outdir + "/evotracer_output/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
#         nwk = outdir + "/evotracer_output/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
#     output:
#         extendedMachina = outdir + "machina_output/" + mouseID + "_all_results_extended.txt",
#         migrationMachina = outdir + "machina_output/" + mouseID + "_migration.txt",
#         seedingMachina = outdir + "machina_output/" + mouseID + "_seeding_topology.txt"
#     params:
#         machinaScripts = "scripts/machina_scripts/",
#         machinaOutPrefix = outdir + "machina_output/" + mouseID,
#         primaryTissue = config['primaryTissue']
#     conda:
#         "envs/machina.yaml"
#     shell:
#         """
#         Module load EBModules;
#         Module load gurobi;

#         scripts/machina_scripts/run_machina.sh --infile {input.asvStat} --tree {input.nwk} --scripts {params.machinaScripts} --prefix {params.machinaOutPrefix} --keep-first-cp
#         """

# rule plotMachinaResults:
#     input:
#         extendedMachina = outdir + "machina_output/" + mouseID + "_all_results_extended.txt",
#         migrationMachina = outdir + "machina_output/" + mouseID + "_migration.txt",
#         seedingMachina = outdir + "machina_output/" + mouseID + "_seeding_topology.txt",
#     output:
#         "test.txt"
#     params:
#         machinaOutDir = outdir + "machina_output"
#     singularity:
#         "envs/evotracer.sif"
#     shell:
#         """
#         Rscript scripts/plotting_scripts/5_machina_analysis/02_machina_tree_v1.R {input.extendedMachina} {params.machinaOutDir};

#         Rscript scripts/plotting_scripts/5_machina_analysis/03_machina_migration_v1.R {input.migrationMachina} {params.machinaOutDir};

#         Rscript scripts/plotting_scripts/5_machina_analysis/04_machina_seeding_topology_v1.R {input.seedingMachina} {params.machinaOutDir}
#         """