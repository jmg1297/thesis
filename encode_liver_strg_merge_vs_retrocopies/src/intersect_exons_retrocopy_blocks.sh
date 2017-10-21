#!/bin/bash

EXONBED=/home/jg600/projects/encode_liver_strg_merge/data/exons/exons.bed

bedtools intersect -wao -a $EXONBED -b $HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.blocks.bed > data/exons_int_retrocopy_blocks/exons.retrocopy_blocks.intersect
