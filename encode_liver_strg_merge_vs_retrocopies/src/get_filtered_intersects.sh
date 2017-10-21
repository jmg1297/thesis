#!/bin/bash

INTERSECTFNAME=data/exons_int_retrocopy_blocks/exons.retrocopy_blocks.intersect
EXONDB=$HOME/projects/encode_liver_strg_merge/data/stringtie_merge/ENCFF001RT_merge.db
RCDB=/home/jg600/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.db
TARGET=data/exons_int_retrocopy_blocks/transcripts.retrocopies.intersect

/usr/bin/python src/filter_block_intersect.py $INTERSECTFNAME $EXONDB $RCDB 0.5 $TARGET

cat $TARGET | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > data/exons_int_retrocopy_blocks/transcripts.retrocopy.bed

cat $TARGET | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > data/exons_int_retrocopy_blocks/expressed_retrocopies.bed

