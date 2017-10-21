#!/bin/bash

SRCDIR=$HOME/projects/cast_rrna_filter/data/rrna_filter
TARGETDIR=data/stringtie

for d in $(find $SRCDIR -type d -name "*CAST*");
do
  SAMPLE=$(echo $d | cut -d/ -f8)
  echo "mkdir "$TARGETDIR"/"$SAMPLE"; stringtie "$d"/rseqc_rrna.ex.sorted.bam -f 0.05 -o "$TARGETDIR"/"$SAMPLE"/out.gtf -M 0.99 -p 4 --fr" >> src/stringtie.cmds
done 
