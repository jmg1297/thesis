#!/bin/bash


SRCDIR=data/strg_merge_int_rt_indels
DBDIR=/home/jg600/projects/stringtie_merge/data/stringtie_merge
RTIBED=/home/jg600/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed

for d in $(find $SRCDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  echo $MRG
  /usr/bin/python src/filter_rti_intersect.py $SRCDIR/$MRG/exons.rti.intersect $DBDIR/$MRG"_merge.db" $RTIBED 0.3 $SRCDIR/$MRG/transcripts.rt_indels.intersect
done
