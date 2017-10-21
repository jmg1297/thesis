#!/bin/bash

SCRIPTDIR=/home/jg600/thesis_projects/strg_merge_all_exons_vs_rt_blocks/src
EXONS=/home/jg600/thesis_projects/encode_liver_strg_merge/data/exons/exons.bed
DATADIR=data/exons_int_rt_blocks

python $SCRIPTDIR/hmap_rt_content.py $DATADIR/transcripts_rt_overlap_summ.bed $EXONS imgs/rt_intersect_hmap.svg 
python $SCRIPTDIR/hmap_rt_content.py $DATADIR/transcripts_rt_overlap_summ.filter_50pc.bed $EXONS imgs/rt_intersect_filt50pc_hmap.svg 
python $SCRIPTDIR/hmap_rt_content.py $DATADIR/transcripts_rt_overlap_summ.filter_90pc.bed $EXONS imgs/rt_intersect_filt90pc_hmap.svg 

