#!/bin/bash

EXONDIR=/home/jg600/projects/stringtie_merge/data/all_exons
RTBLOCKS=/home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed

for d in $(find $EXONDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f8)
  mkdir data/all_exons_int_rt_blocks/$MRG
  bedtools intersect -wao -a $d/exons.bed -b $RTBLOCKS > data/all_exons_int_rt_blocks/$MRG/exons.rt_blocks.intersect
done
