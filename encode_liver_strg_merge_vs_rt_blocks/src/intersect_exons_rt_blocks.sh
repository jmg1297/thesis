#!/bin/bash

EXONS=/home/jg600/thesis_projects/encode_liver_strg_merge/data/exons/exons.bed
RTBLOCKS=/home/jg600/thesis_projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed

mkdir data/exons_int_rt_blocks
bedtools intersect -wao -a $EXONS -b $RTBLOCKS > data/exons_int_rt_blocks/exons.rt_blocks.intersect
