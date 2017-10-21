#!/bin/bash

SCRIPTDIR=/home/jg600/thesis_projects/strg_merge_all_exons_vs_rt_blocks/src
EXONS=/home/jg600/thesis_projects/encode_liver_strg_merge/data/exons/exons.bed

python $SCRIPTDIR/filter_summarised_rt_intersects.py data/exons_int_rt_blocks/transcripts_rt_overlap_summ.bed $EXONS 0.5 data/exons_int_rt_blocks/transcripts_rt_overlap_summ.filter_50pc.bed
python $SCRIPTDIR/filter_summarised_rt_intersects.py data/exons_int_rt_blocks/transcripts_rt_overlap_summ.bed $EXONS 0.9 data/exons_int_rt_blocks/transcripts_rt_overlap_summ.filter_90pc.bed
