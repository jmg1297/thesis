#!/bin/bash

FEATBED=/home/jg600/projects/mm10_ref/rmskJoinedBaseline_blocks+ucscsRetroInfo6_blocks.formatted.bed
ALLRTIBED=/home/jg600/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed
TARGETDIR=/home/jg600/projects/mm10_ref/repeatmasker/retrotransposon_indels/

bedtools intersect -wao \
-a $ALLRTIBED \
-b $FEATBED \
> $TARGETDIR/rt_indels-rmsk+retro_blocks.intersect

python $HOME/projects/src/classify_rt_indels.py \
$TARGETDIR/rt_indels-rmsk+retro_blocks.intersect \
$TARGETDIR/rt_indels.classified.tmp \
30

sort -k1,1 -k2,2n $TARGETDIR/rt_indels.classified.tmp > $TARGETDIR/rt_indels.classified
rm $TARGETDIR/rt_indels.classified.tmp
