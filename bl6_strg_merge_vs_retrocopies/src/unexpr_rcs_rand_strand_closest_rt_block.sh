#!/bin/bash

RMSKFILE=$HOME/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed
BASEDIR=$HOME/projects/strg_merge_all_exons_vs_retrocopies/

for d in $(find data/unexpr_retrocopies_with_rand_strand -type d | sed '1d');
do
echo $d
cd $d
for f in $(find . -type f -name "*.bed");
do
echo $f
bedtools closest -D a -id -io -a $f -b $RMSKFILE > $f.rt_blocks.closest
done
cd $BASEDIR
done
