#!/bin/bash

#####################################
#                                   #
#     Written by Kai Battenberg     #
#     Plant Sciences UC Davis       #
#                                   #
#####################################

#####input options#####
query="Step-02_QUERY/01_QUERY_QCD/QUERY_batch4.fas"
q_seq_type="DNA"
database="Step-04_DATABASE/02_DATABASE/"
db_seq_type="DNA_and_AA"
spp_list="INPUT_NFC/species_list_NFC.txt"
og="Vitis_vinifera"

blast_type="NCBI" #type of BLAST.
w_threshold=3 #minimum word length required for a hit(W in AB-BLAST).
e_threshold="1e-3" #minimum e-value.
id_threshold=0 #minimum %identity.
l_threshold=0 #minimum alignment length.
sander_schneider="YES" #"YES" or "NO", enables threshold described in Sander and Schneider (1991).
lo_threshold=10 #maximum number of loci to keep per species.

i_threshold=1.2 #inflation rate in MCL (between 1.2-6), the higher the I the more things are removed.

rooting="MI" #rooting either as "MD" or "MI".

vraxml="SSE3" #the version of RAxML, AVX or SSE3.
b_threshold=2 #the longest branch length allowed within the last ML tree.

overlap_threshold=0

conversion_list=""

t_threshold=5 #number of threads to use for the run.

toolpath="/share/pilot-kbattenb/bin/:/share/pilot-kbattenb/bin/swipe-2.0.12/Linux/" #path to the tools.
##########



#####Loading modules#####
#module load ab_blast
#module load blast/2.4.0+
#module load swipe
#module load mcl
#module load mafft/7.245
#module load raxml
#module load nw_util
##########



#####Running Orthologfinder#####
#setting up the environment
mkdir -p Step-05_INPUT
mkdir -p Step-05_INPUT/DATABASE
rsync -a ${database}* Step-05_INPUT/DATABASE

rm -f list.txt
list1=$(perl OrthoReD_vNFC/step-05-06_SLURMwrap_step1.pl ${query} ${database})
IFS=',' read -r -a list2 <<< "${list1}"
length=${#list2[@]}
touch list.txt
for i in "${list2[@]}"
do
	echo ${i} >> list.txt
done

#running orthologfinder for each query
qfolder="Step-05_INPUT/individualqueries"
sbatch --array=1-${length} --cpus-per-task=${t_threshold} OrthoReD_vNFC/step-05-06_SLURMwrap_step2.sh\
 list.txt\
 ${q_seq_type}\
 ${database}\
 ${db_seq_type}\
 ${spp_list}\
 ${og}\
 ${blast_type}\
 ${w_threshold}\
 ${e_threshold}\
 ${id_threshold}\
 ${l_threshold}\
 ${sander_schneider}\
 ${lo_threshold}\
 ${i_threshold}\
 ${rooting}\
 ${vraxml}\
 ${b_threshold}\
 ${overlap_threshold}\
 ${t_threshold}\
 ${toolpath}
jobid1=$(squeue -u kbattenb | sort -rnk1 | head -n1 | awk '{print $1}')
IFS='_' read -r -a jobid2 <<< "${jobid1}"

#organizing files
sbatch -c ${t_threshold} -N 1 --mem-per-cpu=2000M --dependency afterok:${jobid2[0]} OrthoReD_vNFC/step-05-06_SLURMwrap_step3.sh\
 ${query}\
 ${q_seq_type}\
 ${database}\
 ${db_seq_type}\
 ${spp_list}\
 ${og}\
 ${blast_type}\
 ${w_threshold}\
 ${e_threshold}\
 ${id_threshold}\
 ${l_threshold}\
 ${sander_schneider}\
 ${lo_threshold}\
 ${i_threshold}\
 ${rooting}\
 ${vraxml}\
 ${b_threshold}\
 ${overlap_threshold}\
 ${t_threshold}\
 ${toolpath}
