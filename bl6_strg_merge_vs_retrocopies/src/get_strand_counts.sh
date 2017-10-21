#!/bin/bash

BASEDIR=$(pwd)

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
cd $d
python $BASEDIR/src/get_strand_counts.py transcript_retrocopy_parent.mbed strand_combo_counts.json
cd $HOME/projects/strg_merge_all_exons_vs_retrocopies
done
