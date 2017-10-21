#!/bin/bash

BASEDIR=$HOME/projects/strg_merge_all_exons_vs_retrocopies/

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
cd $d
comm -23 <(sort /home/jg600/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed) <(sort expressed_retrocopies.bed) | \
sort -k1,1 -k2,2n \
> unexpressed_retrocopies.bed
cd $BASEDIR
done
