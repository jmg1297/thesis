#!/bin/bash

for d in $(find data/novel_int_retrotransposons -type d | sed '1d');
do
samples=$(echo $d | cut -d/ -f3)
echo "python /data/homes/jg600/projects/stringtie_novel_vs_retrotransposons/src/transcript_intersect_hmap.py --intersectFilename "$d"/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect --outputImg imgs/"$samples"_rt_content_hmap.png --sortBy cluster" >> src/rt_content_hmap.cmds
done
