#!/bin/bash

TARGETDIR=data/strg_merge_int_rt_indels
EXONDIR=$HOME/projects/stringtie_merge/data/all_exons
RTIBED=$HOME/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed

for d in $(find $EXONDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f8)
  mkdir $TARGETDIR/$MRG
  bedtools intersect -wao -a $d/exons.bed -b $RTIBED > $TARGETDIR/$MRG/exons.rti.intersect
done
