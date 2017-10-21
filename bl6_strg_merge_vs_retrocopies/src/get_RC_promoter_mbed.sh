#!/bin/bash

TSSMBED=data/ALL_rc_transcript_tss/tss_retrocopy_parent.mbed
GENOMESIZE=$HOME/projects/mm10_ref/genome/mm10.chrom.sizes
TARGET=data/ALL_rc_transcript_tss/promoter_retrocopy_parent.mbed

bedtools slop -l 2000 -r 0 -s -i $TSSMBED -g $GENOMESIZE > $TARGET
