#!/bin/bash

DIR=data/rrna_filter

for d in $(find $DIR -type d -name "*CAST*");
do
  echo "cd "$d"; samtools sort -@ 4 -T tmp.sorted -O 'bam' rseqc_rrna.ex.bam > rseqc_rrna.ex.sorted.bam" >> src/sort_excl_bam.cmds
done
