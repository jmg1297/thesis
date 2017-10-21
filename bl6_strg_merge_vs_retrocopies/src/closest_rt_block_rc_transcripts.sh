#!/bin/bash

TARGETDIR=data/expr_retrocopies_vs_retrotransposon_blocks
RMSKFILE=$HOME/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
MRG=$(echo $d | cut -d/ -f3)
echo $MRG
mkdir $TARGETDIR/$MRG
bedtools closest -D a -id -io -a $d/transcripts.retrocopy.bed -b $RMSKFILE > $TARGETDIR/$MRG/retrocopy_transcripts-rt_blocks.closest
done
