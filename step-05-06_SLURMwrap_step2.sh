#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16000M #RAM per CPU 
#SBATCH --job-name=Ortho_ReD #A single job name for the array
#SBATCH --output=Ortho_ReD_job%A_%a.out #Standard output. "%A" is replaced by the job ID and "%a" with the array index.
#SBATCH --error=Ortho_ReD_job%A_%a.err #Standard error.
#SBATCH --time=8-00:00:00 #maximum run time.

export PATH=${20}:$PATH
S=$(head -n${SLURM_ARRAY_TASK_ID} ${1} | tail -1)

perl OrthoReD_vNFC/step-05-06.pl\
 --query ${S}\
 --q_seq_type ${2}\
 --database ${3}\
 --db_type ${4}\
 --spp_list ${5}\
 --og ${6}\
 --blast_type ${7}\
 --w_threshold ${8}\
 --eval_threshold ${9}\
 --identity_threshold ${10}\
 --length_threshold ${11}\
 --sander_schneider ${12}\
 --loci_threshold ${13}\
 --i ${14}\
 --rooting ${15}\
 --vraxml ${16}\
 --branch_threshold ${17}\
 --overlap_threshold ${18}\
 --procedures 01_02_03_04_05_06_07\
 --threads ${19}\
 --wrap YES
