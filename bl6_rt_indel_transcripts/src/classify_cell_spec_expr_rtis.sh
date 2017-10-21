#!/bin/bash

SCRIPTDIR=$HOME/projects/src
FEATBED=$HOME/projects/mm10_ref/rmskJoinedBaseline_blocks+ucscsRetroInfo6_blocks.formatted.bed
TARGETDIR=data/expr_rti_contents
SRCDIR=data/strg_merge_int_rt_indels

for f in $(find $SRCDIR -type f -name "expr_rtis.*");
do
  name=$(echo $f | cut -d/ -f3 | cut -d. -f1,2)
  echo $f
  echo $name
  bedtools intersect -wao -a $f -b $FEATBED > $TARGETDIR/$name"-rmsk+retro_blocks.intersect"

  python $SCRIPTDIR/classify_rt_indels.py \
  $TARGETDIR/$name"-rmsk+retro_blocks.intersect" \
  $TARGETDIR/$name".classified.tmp" \
  30

  sort -k1,1 -k2,2n $TARGETDIR/$name".classified.tmp" > $TARGETDIR/$name".classified"
  rm $TARGETDIR/$name".classified.tmp"
done
