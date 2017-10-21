#!/bin/bash

for d in $(find data/novel_beds -type d | sed '1d');
do
samples=$(echo $d | cut -d/ -f3)
mkdir data/novel_int_retrotransposons/$samples
echo "bedtools intersect -wo -a "$d"/novel_non-coding.bed -b /home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.formatted.bed > data/novel_int_retrotransposons/"$samples"/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect" >> src/novel_int_rts.cmds
done
