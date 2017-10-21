#!/bin/bash

EXONS=$HOME/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.exons.bed
RTBLOCKS=$HOME/thesis_projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed
TARGET=data/ensembl_exons-rt_blocks.intersect

bedtools intersect -wao -a $EXONS -b $RTBLOCKS > $TARGET
