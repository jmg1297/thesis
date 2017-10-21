#!/bin/bash

for d in $(find data/all_exons_int_retrocopy_blocks/ -type d | sed '1d');
do
/usr/bin/python src/append_parent_bed.py $d/transcripts.retrocopies.filtered-intersect /home/jg600/projects/mm10_ref/retrogenes/retrogene_parent.db $d/transcript_retrocopy_parent.mbed
done
