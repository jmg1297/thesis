#!/bin/bash

RCPROMOTERS=data/ALL_rc_transcript_tss/promoter_retrocopy_parent.mbed
ENSPROMOTERS=$HOME/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.promoters.bed
TARGET=data/ALL_rc_transcript_tss/promoter_retrocopy_parent-ensembl_promoter.intersect

bedtools intersect -wo -a $RCPROMOTERS -b $ENSPROMOTERS > $TARGET
