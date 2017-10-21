#!/bin/bash

SRCDIR=data/all_exons_int_rt_blocks
EXONDIR=/home/jg600/projects/stringtie_merge/data/all_exons
EXONFNAME=exons.bed
RTBLOCKBED=/home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed

for d in $(find $SRCDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  echo $MRG
  python src/collapse_exon_block_intersect.py $d/exons.rt_blocks.intersect $EXONDIR/$MRG/$EXONFNAME $RTBLOCKBED $d/transcripts_retrotransposons.mbed
  python src/summarise_rt_overlaps.py $d/transcripts_retrotransposons.mbed $d/transcripts_rt_overlap_summ.bed
  python src/filter_summarised_rt_intersects.py $d/transcripts_rt_overlap_summ.bed $EXONDIR/$MRG/$EXONFNAME 0.5 $d/transcripts_rt_overlap_summ.filter_50pc.bed
  python src/filter_summarised_rt_intersects.py $d/transcripts_rt_overlap_summ.bed $EXONDIR/$MRG/$EXONFNAME 0.9 $d/transcripts_rt_overlap_summ.filter_90pc.bed
done
