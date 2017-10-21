#!/bin/bash

EXONDIR=/home/jg600/projects/stringtie_merge/data/all_exons

for d in $(find $EXONDIR -type d | sed '1d');
do
MRG=$(echo $d | cut -d/ -f8)
mkdir data/all_exons_int_retrocopy_blocks/$MRG
bedtools intersect -wao -a $d/exons.bed -b $HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.blocks.bed > data/all_exons_int_retrocopy_blocks/$MRG/exons.retrocopy_blocks.intersect
done
