#!/bin/bash

SCRIPTDIR=/home/jg600/thesis_projects/strg_merge_all_exons_vs_rt_blocks/src

python $SCRIPTDIR/summarise_rt_overlaps.py data/exons_int_rt_blocks/transcripts_retrotransposons.mbed data/exons_int_rt_blocks/transcripts_rt_overlap_summ.bed

