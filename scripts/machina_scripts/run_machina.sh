#!/usr/bin/env bash

absolute_path () {
 abspath=`echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"`
 echo $abspath
}


if [[ $# -eq 0 ]] ; then
    echo "Usage: run_machina.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts> --prefix <output_prefix> --primary-tissue <tissue> [--keep-first-cp --threads <int> --batches <int>]"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--infile) ASV="$2"; shift ;;
        -t|--tree) TREE="$2"; shift ;;
        -s|--scripts) SPATH="$2"; shift ;;
        -p|--prefix) PREFIX="$2"; shift ;;
        -o|--primary-tissue) PTISSUE="$2"; shift ;;
        -c|--cutoff) CUTOFF="$2"; shift ;;
        -l|--threads) THREADS=$2; shift ;;
        -b|--batches) BATCHES=$2; shift ;;
        -k|--keep-first-cp) KEEP=true;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: run_machina.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts> --primary-tissue <tissue> [--keep-first-cp --threads <int> --batches <int>]" ; exit 1 ;;
    esac
    shift
done

if ! command -v pmh_tr &> /dev/null
then
    echo "pmh_tr from the MACHINA installation could not be found. Exiting!"
    exit
fi

if [ -z "${CUTOFF}" ]
then
    CUTOFF=1
fi
if [ -z "${THREADS}" ]
then
    THREADS=1
fi
if [ -z "${BATCHES}" ]
then
    BATCHES=1
fi

if [ -z "${ASV}" ]
then
    echo "Missing --infile parameter. Exiting!"
    exit
fi
if [ -z "${TREE}" ]
then
     echo "Missing --tree parameter. Exiting!"
    exit
fi
if [ -z "${SPATH}" ]
then
    echo "Missing --scripts parameter. Exiting!"
    exit
fi
if [ -z "${PREFIX}" ]
then
    echo "Missing --prefix parameter. Exiting!"
    exit
fi
if [ -z "${PTISSUE}" ]
then
    echo "Missing --primary-tissue parameter. Exiting!"
    exit
fi

# PATHS TO SCRIPTS 
BIG_CP_THRESHOLD="${CUTOFF//[$'\t\r\n ']}"

SPATH=`absolute_path "${SPATH}"` 
TREE=`absolute_path "${TREE}"`
ASV=`absolute_path "${ASV}"`
TRAV="${SPATH}/traverse_split.py"
GET="${SPATH}/get_results.sh"
GETOLD="${SPATH}/get_results_old2new.sh"
MACHINA="pmh_tr"
TOPOLOGY="${SPATH}/print_seeding_topology.py"
SELECTION="${SPATH}/selection_tree_test.py"
MIGRATION="${SPATH}/count_migrations.py"
ADD_INFO="${SPATH}/add_freq_prob_to_results.py"
# Ensure this harcoded reference
REFERENCE="BC10"

## PREPROCESS INPUT DATA ##

if [[ -n "${BIG_CP_THRESHOLD//[0-9]/}" ]]; then
    echo "Value for cutoff parameter is not an integer!"
    exit
fi

if [ -d "${PREFIX}/cp_output" ]; then
  echo "Output directory ${PREFIX}/cp_output already exist. Exiting!"
  exit
fi

mkdir ${PREFIX}/cp_output
cd ${PREFIX}/cp_output
# Extract key ASV columns (asv_names,sample,group)
cut -d',' -f1,2,30 ${ASV} | grep -v "$REFERENCE"| sed '/^$/d' > ${PREFIX}/asv_sample_group.csv
# exclude miscelleneaous group CP00
if [ "$KEEP" = true ] ; then
    remove="^$"
    ###remove="^CP01$" ### ONLY FOR TESTING TO REMOVE CP01 FOR SPEED
else
    # remove CP0, which can be zero padded, e.g. "CP00", thus we sort
    first_cp=`cut -f3 -d',' ${PREFIX}/asv_sample_group.csv|tail -n +2| sort| head -1`
    remove="^${first_cp}$"
fi
# if a CP only includes the reference, exclude it
cut -f3 -d',' ${PREFIX}/asv_sample_group.csv|tail -n +2| sort| uniq | egrep -v "$remove"| sed '/^$/d' > ${PREFIX}/CP_list.txt
# if its empty, then exit
if [ ! -s ${PREFIX}/CP_list.txt ]; then
    echo "No CPs found or all CPs removed after filtering. Exiting!"
    exit
fi

# prepare machina input files for each CP
touch ${PREFIX}/big_CP_list.txt
while read l; do
    ${TRAV} ${TREE} ${PREFIX}/asv_sample_group.csv ${l} ${PTISSUE}
    num_labels=`sed -n '$=' ${l}_labels_split.txt`
    if [[ $num_labels -gt $BIG_CP_THRESHOLD ]]; then
    #if [ "$num_labels" -gt 30 ]; then
        echo ${l} >> ${PREFIX}/big_CP_list.txt
    fi
done<${PREFIX}/CP_list.txt
# remove CPs that failed traversal due to unexpected clade topology
#grep -v -f FailedCP.txt ${PREFIX}/CP_list.txt > ${PREFIX}/CP_list.tmp
#mv ${PREFIX}/CP_list.tmp ${PREFIX}_CP/list.txt
# make output directory for each CP
while read l; do mkdir ${l}_split;done<${PREFIX}/CP_list.txt
   
# prepare input for selection scan on original tree
for l in *_tree_split.txt; do cat $l |sed "s/^/${l} tree /"| sed 's/_tree_split.txt//';done | grep -v "$remove" > ${PREFIX}/all_original_tree.txt

## RUN MACHINA ##     
# make command file for machina
# speed up MACHINA a bit by adding "-m 3"
grep -f ${PREFIX}/big_CP_list.txt ${PREFIX}/CP_list.txt| while read l; do echo "${MACHINA} -OLD -t ${THREADS} -m 3 -o ${l}_split -c ${l}_colors.txt -p ${PTISSUE} ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done >> ${PREFIX}/machina.cmd
grep -v -f ${PREFIX}/big_CP_list.txt ${PREFIX}/CP_list.txt| while read l; do echo "${MACHINA} -t ${THREADS} -m 3 -o ${l}_split -c ${l}_colors.txt -p ${PTISSUE} ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done >> ${PREFIX}/machina.cmd

# run machina in parallel
module load EBModules
module load Gurobi
ParaFly -CPU ${BATCHES} -c ${PREFIX}/machina.cmd
module unload EBModules
module unload Gurobi

# parse results from each machina output dir
grep -f ${PREFIX}/big_CP_list.txt ${PREFIX}/CP_list.txt| while read l; do ${GETOLD} $l ${PTISSUE} ${SPATH};done | tr '\t' ' '>> ${PREFIX}/all_results.txt
grep -v -f ${PREFIX}/big_CP_list.txt ${PREFIX}/CP_list.txt| while read l; do ${GET} $l ${PTISSUE} ${SPATH};done | tr '\t' ' '>> ${PREFIX}/all_results.txt

## ANALYSE INFERRED TOPOLOGY
python $TOPOLOGY ${PREFIX}/all_results.txt ${PTISSUE} > ${PREFIX}/seeding_topology.txt 
python $MIGRATION ${PREFIX}/all_results.txt > ${PREFIX}/migration.txt

## ANALYSE SELECTION ON ORIGINAL AND MACHINA TOPOLOGY

#python $SELECTION ${PREFIX}/cp_output/${PREFIX}/all_results.txt $ASV > ${PREFIX}/cp_output/${PREFIX}/selection.txt
python $SELECTION ${PREFIX}/all_original_tree.txt $ASV| grep "^test" > ${PREFIX}/selection_original_test.txt
python $SELECTION ${PREFIX}/all_original_tree.txt $ASV| grep "^expansion" > ${PREFIX}/selection_original_expansion.txt
python $ADD_INFO ${PREFIX}/migration.txt $ASV ${PREFIX}/all_results.txt > ${PREFIX}/all_results_extended.txt

# Clean up
cd ..
mkdir data
mv *cp_output data/
mv *CP_list.txt data/
mv *all_results.txt data/
mv *all_original_tree.txt data/
mv *asv_sample_group.csv data/
mv *.cmd data/
mv *.cmd.completed data/
