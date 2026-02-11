import os, sys
from glob import glob

# get mouse ids and fastq files from input fastq dir
mice = []
mouse_fastqs = {}
for file in glob(os.path.join(config['fastqPath'], "*.fastq.gz")):
    mouse = os.path.basename(file).split("_")[0]
    if mouse not in mice:
        mice.append(mouse)
    if mouse not in mouse_fastqs:
        mouse_fastqs[mouse] = []
    mouse_fastqs[mouse].append(file)
    

def get_elements_from_file(path):
    # to get all cps
    filename = f"{path}/cp_list.txt"
    with open(filename) as f:
        cps = [line.strip() for line in f if set(list(line.strip().replace("CP", ""))) != set('0')]

    cps = [cp for cp in cps if cp != "CP01"]

    return cps


outdir = config['outdir']
envs = config['envs']



rule all:
    input:
        # expand("{outdir}/evotracer/{mouse}/stat_cps_dispersal_bargraph_hm.pdf", outdir = outdir, mouse = mice),
        expand("{outdir}/cp_split/{mouse}/cp_list.txt", outdir = outdir, mouse = mice),
        expand("{outdir}/cp_split/{mouse}/cp_size_barplot.pdf", outdir = outdir, mouse = mice),

        # the above are required to run first with the below files commented out to first generate the cp list which can then be called by those below
        [expand("{outdir}/mach2/{mouse}/{cp}/consensus_graph_filtered_by_threshold.pdf", outdir = outdir, mouse = m, cp = get_elements_from_file(f"{outdir}/cp_split/{m}")) for m in mice],
        [expand("{outdir}/mach2/{mouse}/{cp}/all_transition_matrix.pdf", outdir = outdir, mouse = m, cp = get_elements_from_file(f"{outdir}/cp_split/{m}")) for m in mice],
        [expand("{outdir}/mach2/{mouse}/{cp}/all_seeding_topologies.pdf", outdir = outdir, mouse = m, cp = get_elements_from_file(f"{outdir}/cp_split/{m}")) for m in mice],
        expand("{outdir}/mach2/{mouse}/overall_transition_matrix.pdf", outdir=outdir, mouse=mice),
        expand("{outdir}/mach2/{mouse}/overall_seeding_topologies.pdf", outdir=outdir, mouse=mice),

rule runEvotracer:
    input:
        fastqFiles = lambda wildcards: mouse_fastqs[wildcards.mouse]
    output:
        char_matrix = "{outdir}/evotracer/{mouse}/asv_analysis/character_matrix.csv",
        asv_stat = "{outdir}/evotracer/{mouse}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
        mut_dict = "{outdir}/evotracer/{mouse}/asv_analysis/mutation_profile_map.csv",
        evotracerCassiopeiaTree= "{outdir}/evotracer/{mouse}/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
        rDataObject = "{outdir}/evotracer/{mouse}/evotracer.RData"
    params:
        asvThreshold = config['asvThreshold'],
        refName = config['ref_name'],
        refSeq = config['ref_seq'],
        refFlankLeft = config['ref_flank_left'],
        refFlankRight = config['ref_flank_right'],
        refCutSites = config['ref_cut_sites'],
        refBorderSites = config['ref_border_sites'],
        evoOutDir = "{outdir}/evotracer/{mouse}",
        scripts = config['scripts'],
        threads = 10,
        mem = '1G',
    singularity:
        f"{envs}/evotracer.sif"
    shell:
        """
        # make mouse specific temp dir with fastqs
        tempDir={wildcards.outdir}/{wildcards.mouse}_fastq_temp
        mkdir -p $tempDir

        for file in {input.fastqFiles}; do
            cp $file $tempDir/
        done

        # evotracer needs fastqs to be unzipped
        gunzip $tempDir/*

        Rscript {params.scripts}/evotracer_scripts/evotracer.R $tempDir \
        "{params.evoOutDir}" \
        "{params.asvThreshold}" \
        "{params.refName}" \
        "{params.refSeq}" \
        "{params.refFlankLeft}" \
        "{params.refFlankRight}" \
        "{params.refCutSites}" \
        "{params.refBorderSites}"

        # remove temp dir
        rm -r $tempDir
        """

rule plotEvotracerResults:
    input:
        rDataObject = "{outdir}/evotracer/{mouse}/evotracer.RData",
    output:
        indels = "{outdir}/evotracer/{mouse}/hist_freq_indels.pdf",
        dispersal = "{outdir}/evotracer/{mouse}/stat_cps_dispersal_bargraph_hm.pdf",
        phylo = "{outdir}/evotracer/{mouse}/cp_tree_msa_cna_bc_bubble_qnt_ggtree_mp.pdf",
    params:
        evoPlotsOutDir = "{outdir}/evotracer/{mouse}/evotracer_plots",
        scripts = config['scripts'],
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/evotracer_plotting.sif"
    shell:
        """
        Rscript {params.scripts}/plotting_scripts/2_barcode_edits_analysis/04.1_hist_freq_indels.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript {params.scripts}/plotting_scripts/2_barcode_edits_analysis/04.2_hist_freq_seq_length.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript {params.scripts}/plotting_scripts/2_barcode_edits_analysis/04.3_hist_freq_site_marked.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript {params.scripts}/plotting_scripts/3_clonal_population_analysis/05.1_cp_stat_vis.R {input.rDataObject} {params.evoPlotsOutDir};

        Rscript {params.scripts}/plotting_scripts/4_phylogenetic_analysis/06.1_tree_msa_bubble_all_clones.R {input.rDataObject} {params.evoPlotsOutDir}
        """

rule splitCPs:
    input:
        char_matrix = "{outdir}/evotracer/{mouse}/asv_analysis/character_matrix.csv",
        asv_stat = "{outdir}/evotracer/{mouse}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv",
        mut_dict = "{outdir}/evotracer/{mouse}/asv_analysis/mutation_profile_map.csv",
    output:
        cpList = "{outdir}/cp_split/{mouse}/cp_list.txt",
    params:
        outdir = "{outdir}/cp_split/{mouse}",
        scripts = config['scripts'],
        threads = 1,
        mem = '1G',
    shell:
        """
        # prep mutation dictionary for the mouse overall
        sed 's/\",/\\t/g' {input.mut_dict} | sed 's/\"//g' > {params.outdir}/mutation_dict.tsv

        cps=$(tail -n +2 {input.asv_stat} | cut -d',' -f30 | sort | uniq | tr '\\n' ' ')
        for cp in $cps; do

            # make cp dir
            cpOutdir={params.outdir}/${{cp}}
            mkdir -p $cpOutdir

            # get ASVs for the cp
            asvs=$(tail -n +2 {input.asv_stat} | grep -w "$cp" | cut -d',' -f1 | sort | uniq | tr '\\n' ' ')

            # setup temp matrix and tissues headers
            tempMatrix=$cpOutdir/matrix.tsv
            head -n 1 {input.char_matrix} | sed 's/,/\\t/g' | sed 's/\"//g' > $tempMatrix
            tempTissues=$cpOutdir/tissues.tsv
            echo -e "group_name\\ttissues" > $tempTissues

            # subset matrix and tissues for each cp to only retain ASVs for that CP
            for asv in $asvs; do
                sed 's/,/\\t/g' {input.char_matrix} | sed 's/\"//g' | grep -w "$asv" >> $tempMatrix
                tail -n +2 {input.asv_stat} | cut -d',' -f1-2 | grep -w "$asv" | awk -F',' '{{if (!seen[$1","$2]++) a[$1]=a[$1]","$2}} END {{for (i in a) {{split(a[i], b, ","); delete seen; c=""; for (j in b) if (b[j] != "" && !seen[b[j]]++) c=(c == "" ? b[j] : c","b[j]); print i"\\t"c}}}}' >> ${{tempTissues}}
            done

            # track all CPs as they finish pre-processing
            echo "$cp" >> {output.cpList}
        done
        """

rule plotCpSizeBars:
    input:
        cpLists = "{outdir}/cp_split/{mouse}/cp_list.txt"
    output:
        asvCountPdf = "{outdir}/cp_split/{mouse}/cp_size_barplot.pdf"
    params:
        outdir = "{outdir}/cp_split/{mouse}",
        scripts = config['scripts'],
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/networkx.sif"
    shell:
        """
        # Count the number of ASVs in each CP for each mouse
        asv_count_csv={params.outdir}/asv_counts_per_cp.csv
        echo "mouse,cp,num_asvs" > $asv_count_csv

        files=$(find {params.outdir} -type f -name "tissues.tsv")

        for file in $files; do
            mouse=$(echo $file | rev | cut -d'/' -f3 | rev)
            cp=$(echo $file | rev | cut -d'/' -f2 | rev)
            num_asvs=$(($(wc -l < $file) - 1))
            echo "$mouse,$cp,$num_asvs" >> $asv_count_csv
        done

        # Plot the number of ASVs in each CP for each mouse
        python {params.scripts}/plotting_scripts/plot_bar_cp_sizes.py $asv_count_csv {output.asvCountPdf}
        """


rule prepMachina:
    input:
        evotracerCassiopeiaTree = "{outdir}/evotracer/{mouse}/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick",
        cpList = "{outdir}/cp_split/{mouse}/cp_list.txt",
    output:
        labeling = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}.labeling",
        colors = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}_colors.txt",
        tree = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}.tree",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        outputdir = "{outdir}/machina_prep/{mouse}/{cp}",
        outprefix = "{mouse}_{cp}",
        outdir = "{outdir}",
        scripts = config['scripts'],
        threads = 1,
        mem = '1G'
    singularity:
        f"{envs}/simulate.sif"
    shell:
        """
        tissuesTsv={wildcards.outdir}/cp_split/{wildcards.mouse}/{wildcards.cp}/tissues.tsv

        # modify the tissues so that there is always only one primaryTissue name, so if the primaryTissue name is found in the tissues with additional A,B or 1,2 etc. labels then will all be turned to primary label only

        # prep inputs
        python {params.scripts}/machina_scripts/prep_machina.py {input.evotracerCassiopeiaTree} {params.primaryTissue} $tissuesTsv {params.outputdir} {params.outprefix}
        """


rule runMach2:
    input:
        labeling = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}.labeling",
        colors = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}_colors.txt",
        tree = "{outdir}/machina_prep/{mouse}/{cp}/{mouse}_{cp}.tree",
    output:
        mach2graph = "{outdir}/mach2/{mouse}/{cp}/BDR-G-0.graph",
    params:
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        threads = 50,
        mem = '1G',
    conda:
        f"{envs}/mach2.yaml"
    shell:
        """
        unique_tissues=$(cut -f2 {input.labeling} | sort | uniq)
        num_unique_tissues=$(echo "$unique_tissues" | wc -l)

        if [ "$num_unique_tissues" -gt 1 ]; then
            # check if primary tissue is in the labeling file
            if grep -q "{params.primaryTissue}" {input.labeling}; then
                mach2 {input.tree} {input.labeling} --colormap {input.colors} -p {params.primaryTissue} --log -o {params.outputdir} -t {params.threads}

                # Need to touch output files to satisfy output, only if output with a different name exists
                if [ $(ls {params.outputdir}/*.graph | wc -l) -gt 0 ] && [ ! -f {output.mach2graph} ]; then
                    touch {output}
                fi
            else
                most_common_tissue=$(cut -f2 {input.labeling} | sort | uniq -c | sort -nr | head -n 1 | awk '{{print $2}}')
                echo "The {params.primaryTissue} is not in the tissue labels. Using the most common tissue $most_common_tissue as the primary tissue."
                mach2 {input.tree} {input.labeling} --colormap {input.colors} -p $most_common_tissue --log -o {params.outputdir} -t {params.threads}

                # Need to touch output files to satisfy output, only if output with a different name exists
                if [ $(ls {params.outputdir}/*.graph | wc -l) -gt 0 ] && [ ! -f {output.mach2graph} ]; then
                    touch {output}
                fi
            fi
        else
            echo "Only one unique tissue found in the labeling file, so not enough to run Machina. Exiting."
            touch {output}
        fi

        if [ $(ls {params.outputdir}/*.graph | wc -l) -gt 0 ] && [ ! -f {output.mach2graph} ]; then
            touch {output}
        fi
        """

rule gatherMach2Results:
    input:
        mach2graph = "{outdir}/mach2/{mouse}/{cp}/BDR-G-0.graph"    # Not sure how to automate the primary tissue name in the input
    output:
        mach2Combined = "{outdir}/mach2/{mouse}/{cp}/graph_results_combined.txt",
        mach2Consensus = "{outdir}/mach2/{mouse}/{cp}/consensus_graph.txt",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    shell:
        """
        # get all graph files
        results=$(find $(dirname {input.mach2graph}) -name "*.graph")

        # get all unique tissue names in results where the file basename is tissue-*
        unique_tissues=$(for file in $results; do basename $file | cut -d'-' -f1; done | sort | uniq)
        useRootTissue=false
        if [ $(echo "$unique_tissues" | wc -l) -gt 1 ]; then
            root_tissue=$(echo "$unique_tissues" | grep -v "{params.primaryTissue}" | head -n 1)
            useRootTissue=true
        fi

        # combine all results
        echo "result_num,source,target,count" > {output.mach2Combined}
        for file in $results; do
            if [ -s "$file" ]; then
                result_num=$(basename $file | cut -d'-' -f3 | cut -d'.' -f1)
                if [ "$useRootTissue" = true ]; then
                    echo "$result_num,{params.primaryTissue},$root_tissue,1" >> {output.mach2Combined}
                fi
                while IFS= read -r line; do
                    echo "$result_num,$(echo $line | tr ' ' ',')" >> {output.mach2Combined}
                done < "$file"
            fi
        done

        # get the probability weighted edge consensus graph
        python {params.scripts}/machina_scripts/convert_mach2_results_to_consensus_graph.py {output.mach2Combined} {output.mach2Consensus}
        """

rule filterMach2Consensus:
    input:
        mach2Consensus = "{outdir}/mach2/{mouse}/{cp}/consensus_graph.txt",
    output:
        mach2Filtered = "{outdir}/mach2/{mouse}/{cp}/consensus_graph_filtered_by_threshold.txt",
    params:
        consensusThreshold = 0.5,
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    shell:
        """
        for line in $(cat {input.mach2Consensus}); do
            freq=$(echo $line | cut -d',' -f2)
            if (( $(echo "$freq > {params.consensusThreshold}" | bc -l) )); then
                echo $line >> {output.mach2Filtered}
            fi
        done

        if [ ! -s {output.mach2Filtered} ]; then
            touch {output.mach2Filtered}
            echo "The output file {output.mach2Filtered} was empty, so no migration graph edges passed the consensus threshold set."
        fi
        """

rule plotFilteredMach2ConsensusGraph:
    input:
        mach2Filtered = "{outdir}/mach2/{mouse}/{cp}/consensus_graph_filtered_by_threshold.txt",
    output:
        mach2Graph = "{outdir}/mach2/{mouse}/{cp}/consensus_graph_filtered_by_threshold.pdf",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    conda:
        f"{envs}/networkx3.yaml"
    shell:
        """
        python {params.scripts}/plotting_scripts/plot_mach2_consensus_graph.py {input.mach2Filtered} {params.primaryTissue} {output.mach2Graph}
        """

# rule plotMach2AllGraphsAndTrees:
#     input:
#         mach2graph = "{outdir}/mach2/{mouse}/{cp}/BDR-G-0.graph"
#     output:
#         mach2graphPlot = "{outdir}/mach2/{mouse}/{cp}/BDR-G-0.pdf"
#     params:
#         primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
#         scripts = config['scripts'],
#         outputdir = "{outdir}/mach2/{mouse}/{cp}",
#         threads = 1,
#         mem = '1G',
#     conda:
#         f"{envs}/networkx3.yaml"
#     shell:
#         """
#         files=$(find {params.outputdir} -name "*.dot")

#         for file in $files; do
#             outfile=$(echo $file | sed 's/.dot/.pdf/')
#             python {params.scripts}/plotting_scripts/plot_migration_graph_from_dot.py $file $outfile
#         done

#         # Need to touch output files to satisfy output, only if output with a different name exists
#         if [ $(ls {params.outputdir}/*-G-*.pdf | wc -l) -gt 0 ] && [ ! -f {output.mach2graphPlot} ]; then
#             touch {output.mach2graphPlot}
#         elif [ $(ls {params.outputdir}/*-G-*.dot | wc -l) -eq 0 ]; then
#             touch {output.mach2graphPlot}
#         fi
#         """

rule callTransitionMatricesPerCP:
    input:
        mach2graph = "{outdir}/mach2/{mouse}/{cp}/BDR-G-0.graph"
    output:
        allSolutionsTransitionMatrix = "{outdir}/mach2/{mouse}/{cp}/all_transition_matrix.csv",
        allSeedingTopologies = "{outdir}/mach2/{mouse}/{cp}/all_seeding_topologies.csv",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    conda:
        f"{envs}/plotting.yaml"
    shell:
        """
        # use every tree file to call a transition matrix
        files=$(find {params.outputdir} -name "*.tree")

        if [ -z "$files" ]; then
            echo "No solutions exist, so not plotting transition matrices."
            touch {output}
        else
            for file in $files; do
                if [[ "$file" == *".refined.tree" ]]; then
                    labelingFile=$(echo $file | sed 's/.refined.tree/.location.labeling/')
                    outfile=$(echo $file | sed 's/.refined.tree/_transition_matrix.csv/')
                    outfile2=$(echo $file | sed 's/.refined.tree/_seeding_topologies.csv/')
                else
                    labelingFile=$(echo $file | sed 's/.tree/.labeling/')
                    outfile=$(echo $file | sed 's/.tree/_transition_matrix.csv/')
                    outfile2=$(echo $file | sed 's/.tree/_seeding_topologies.csv/')
                fi
                python {params.scripts}/machina_scripts/call_transition_matrix_from_mach2_result.py $file $labelingFile $outfile {params.primaryTissue}
                python {params.scripts}/machina_scripts/get_seeding_topologies.py $file $labelingFile {params.primaryTissue} $outfile2 {wildcards.cp}
            done

            transitionMatrixFiles=$(find {params.outputdir} -name "*_transition_matrix.csv" | paste -sd "," -)

            if [ -n "$transitionMatrixFiles" ]; then
                python {params.scripts}/machina_scripts/combine_transition_matrices.py "$transitionMatrixFiles" "{output.allSolutionsTransitionMatrix}" "False"
            fi


            seedingTopologyFiles=$(find {params.outputdir} -name "*_seeding_topologies.csv" | paste -sd "," -)

            if [ -n "$seedingTopologyFiles" ]; then
                python {params.scripts}/machina_scripts/combine_seeding_topologies.py "$seedingTopologyFiles" "{output.allSeedingTopologies}" "False"
            fi
        fi
        """

rule plotTransitionMatricesPerCP:
    input:
        allSolutionsTransitionMatrix = "{outdir}/mach2/{mouse}/{cp}/all_transition_matrix.csv",
    output:
        allSolutionsTransitionMatrixPlot = "{outdir}/mach2/{mouse}/{cp}/all_transition_matrix.pdf",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/evotracer_plotting.sif"
    shell:
        """
        transitionMatrixFiles=$(find {params.outputdir} -name "*_transition_matrix.csv")

        for file in $transitionMatrixFiles; do
        outfile=$(echo $file | sed 's/.csv/.pdf/')
            if [ -s "$file" ]; then
                Rscript {params.scripts}/plotting_scripts/5_machina_analysis/03_machina_migration_v1.R "$file" "$outfile" "{params.primaryTissue}"
            else
                touch $outfile
            fi
        done
        """

rule plotSeedingTopologyPerCP:
    input:
        allSeedingTopologies = "{outdir}/mach2/{mouse}/{cp}/all_seeding_topologies.csv",
    output:
        allSeedingTopologiesPlot = "{outdir}/mach2/{mouse}/{cp}/all_seeding_topologies.pdf",
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}/{cp}",
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/evotracer_plotting.sif"
    shell:
        """
        seedingTopologyFiles=$(find {params.outputdir} -name "*_seeding_topologies.csv")

        for file in $seedingTopologyFiles; do
            outfile=$(echo $file | sed 's/.csv/.pdf/')
            if [ -s "$file" ]; then
                Rscript {params.scripts}/plotting_scripts/5_machina_analysis/04_machina_seeding_topology_v1.R "$file" "$outfile"
            else
                touch $outfile
            fi
        done
        """


rule getOverallTransitionMatrixPerMouse:
    input:
        lambda wildcards: expand("{outdir}/mach2/{mouse}/{cp}/all_transition_matrix.csv", outdir='{outdir}', mouse='{mouse}', cp=get_elements_from_file(f"{wildcards.outdir}/cp_split/{wildcards.mouse}")),
        # "{outdir}/mach2/{mouse}/CP03/all_transition_matrix.csv"
    output:
        overallTransitionMatrix = "{outdir}/mach2/{mouse}/overall_transition_matrix.csv",
    params:
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}",
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/networkx.sif"
    shell:
        """
        transitionMatrixFiles=$(echo "{input}" | sed 's/ /,/g')

        python {params.scripts}/machina_scripts/combine_transition_matrices.py "$transitionMatrixFiles" {output.overallTransitionMatrix} True
        """

rule plotOverallTransitionMatrixPerMouse:
    input:
        overallTransitionMatrix = "{outdir}/mach2/{mouse}/overall_transition_matrix.csv",
    output:
        overallTransitionMatrixPlot = "{outdir}/mach2/{mouse}/overall_transition_matrix.pdf"
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}",
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/evotracer_plotting.sif"
    shell:
        """
        Rscript {params.scripts}/plotting_scripts/5_machina_analysis/03_machina_migration_v1.R "{input.overallTransitionMatrix}" "{output.overallTransitionMatrixPlot}" "{params.primaryTissue}"
        """

rule getOverallSeedingTopologyPerMouse:
    input:
        lambda wildcards: expand("{outdir}/mach2/{mouse}/{cp}/all_seeding_topologies.csv", outdir='{outdir}', mouse='{mouse}', cp=get_elements_from_file(f"{wildcards.outdir}/cp_split/{wildcards.mouse}")),
        # "{outdir}/mach2/{mouse}/CP03/all_seeding_topologies.csv"
    output:
        overallSeedingTopologies = "{outdir}/mach2/{mouse}/overall_seeding_topologies.csv",
    params:
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}",
        threads = 1,
        mem = '1G',
    conda:
        f"{envs}/plotting.yaml"
    shell:
        """
        seedingTopologyFiles=$(echo "{input}" | sed 's/ /,/g')

        python {params.scripts}/machina_scripts/combine_seeding_topologies.py "$seedingTopologyFiles" "{output.overallSeedingTopologies}" "True"
        """

rule plotOverallSeedingTopologyPerMouse:
    input:
        overallSeedingTopologies = "{outdir}/mach2/{mouse}/overall_seeding_topologies.csv",
    output:
        overallSeedingTopologiesPlot = "{outdir}/mach2/{mouse}/overall_seeding_topologies.pdf"
    params:
        primaryTissue = lambda wildcards: config['primaryTissue'][wildcards.mouse],
        scripts = config['scripts'],
        outputdir = "{outdir}/mach2/{mouse}",
        threads = 1,
        mem = '1G',
    singularity:
        f"{envs}/evotracer_plotting.sif"
    shell:
        """
        Rscript {params.scripts}/plotting_scripts/5_machina_analysis/04_machina_seeding_topology_v1.R "{input.overallSeedingTopologies}" "{output.overallSeedingTopologiesPlot}"
        """

