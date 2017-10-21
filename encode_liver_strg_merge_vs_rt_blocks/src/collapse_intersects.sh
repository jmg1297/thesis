#!/bin/bash

SCRIPTDIR=/home/jg600/thesis_projects/strg_merge_all_exons_vs_rt_blocks/src
INTERSECT=data/exons_int_rt_blocks/exons.rt_blocks.intersect
EXONS=/home/jg600/thesis_projects/encode_liver_strg_merge/data/exons/exons.bed
RTBLOCKS=/home/jg600/thesis_projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed


python $SCRIPTDIR/collapse_exon_block_intersect.py $INTERSECT $EXONS $RTBLOCKS data/exons_int_rt_blocks/transcripts_retrotransposons.mbed

