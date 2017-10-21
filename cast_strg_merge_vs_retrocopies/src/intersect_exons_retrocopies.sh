#!/bin/bash

EXONDIR=/home/jg600/projects/cast_strg_merge/data/exons
RCBED=/home/jg600/projects/cast_retrocopy_alignment/data/filter_by_chrom_length_rel_pos/conserved_retrocopies.blocks.bed
TARGETDIR=data/exons_int_retrocopies

for d in $(find $EXONDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f8)
  mkdir $TARGETDIR/$MRG
  bedtools intersect -wao -a $d/exons.bed -b $RCBED > $TARGETDIR/$MRG/exons.retrocopies.intersect
done
