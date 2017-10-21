#!/bin/bash

SRCDIR=data/all_exons_int_rt_blocks
EXONDIR=/home/jg600/projects/stringtie_merge/data/all_exons
EXONFNAME=exons.bed

for d in $(find $SRCDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  echo "python src/hmap_rt_content.py "$d"/transcripts_rt_overlap_summ.bed "$EXONDIR"/"$MRG"/"$EXONFNAME" imgs/"$MRG"_rt_intersect_hmap.svg" >> src/hmap.cmds
  echo "python src/hmap_rt_content.py "$d"/transcripts_rt_overlap_summ.filter_50pc.bed "$EXONDIR"/"$MRG"/"$EXONFNAME" imgs/"$MRG"_rt_intersect_filt50pc_hmap.svg" >> src/hmap.cmds
  echo "python src/hmap_rt_content.py "$d"/transcripts_rt_overlap_summ.filter_90pc.bed "$EXONDIR"/"$MRG"/"$EXONFNAME" imgs/"$MRG"_rt_intersect_filt90pc_hmap.svg" >> src/hmap.cmds
done
