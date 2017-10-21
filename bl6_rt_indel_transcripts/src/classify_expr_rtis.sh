#!/bin/bash

SCRIPTDIR=$HOME/projects/src
FEATBED=$HOME/projects/mm10_ref/rmskJoinedBaseline_blocks+ucscsRetroInfo6_blocks.formatted.bed
SRCDIR=data/strg_merge_int_rt_indels
TARGETDIR=data/expr_rti_contents

for d in $(find $SRCDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  mkdir $TARGETDIR/$MRG

  bedtools intersect -wao \
  -a $d/expressed_rtis.bed \
  -b $FEATBED \
  > $TARGETDIR/$MRG/expressed_rtis-rmsk+retro_blocks.intersect

  python $SCRIPTDIR/classify_rt_indels.py \
  $TARGETDIR/$MRG/expressed_rtis-rmsk+retro_blocks.intersect \
  $TARGETDIR/$MRG/expressed_rtis.classified.tmp \
  30

  sort -k1,1 -k2,2n $TARGETDIR/$MRG/expressed_rtis.classified.tmp > $TARGETDIR/$MRG/expressed_rtis.classified
  rm $TARGETDIR/$MRG/expressed_rtis.classified.tmp
done
