#!/bin/bash

#SBATCH --output=Ortho_ReD_job%A.out #Standard output. "%A" is replaced by the job ID.
#SBATCH --time=8-00:00:00 #maximum run time.

export PATH=${20}:$PATH

perl OrthoReD_vNFC/step-05-06.pl\
 --query ${1}\
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
 --procedures 08_09\
 --threads ${19}\
 --wrap YES

rm -f list.txt
