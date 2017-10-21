#!/bin/bash

INTERSECT=data/ensembl_exons-rt_blocks.intersect
EXONS=$HOME/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.exons.bed
RTBLOCKS=$HOME/thesis_projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed
TARGET=data/ensembl_transcripts-retrotransposons.mbed

python src/collapse_exon_block_intersect.py $INTERSECT $EXONS $RTBLOCKS $TARGET
