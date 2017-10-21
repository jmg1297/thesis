#!/bin/bash

for d in $(find data/all_exons_int_retrocopy_blocks/ -type d | sed '1d');
do
MRG=$(echo $d | cut -d/ -f3)
/usr/bin/python src/filter_block_intersect.py $d/exons.retrocopy_blocks.intersect /home/jg600/projects/stringtie_merge/data/stringtie_merge/$MRG"_merge.db" /home/jg600/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.db 0.5 $d/transcripts.retrocopies.filtered-intersect

cat $d/transcripts.retrocopies.filtered-intersect | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $d/transcripts.retrocopy.bed

cat $d/transcripts.retrocopies.filtered-intersect | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $d/expressed_retrocopies.bed

done
