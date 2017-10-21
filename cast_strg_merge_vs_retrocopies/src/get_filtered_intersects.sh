#!/bin/bash

SCRIPTDIR=$HOME/projects/strg_merge_all_exons_vs_retrocopies/src
DBDIR=$HOME/projects/cast_strg_merge/data/stringtie_merge
RCDB=$HOME/projects/cast_retrocopy_alignment/data/filter_by_chrom_length_rel_pos/conserved_retrocopies.db

for d in $(find data/exons_int_retrocopies -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  /usr/bin/python $SCRIPTDIR/filter_block_intersect.py $d/exons.retrocopies.intersect $DBDIR/$MRG"_merge.db" $RCDB 0.5 $d/transcripts.retrocopies.filtered-intersect

  cat $d/transcripts.retrocopies.filtered-intersect | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $d/transcripts.retrocopy.bed

  cat $d/transcripts.retrocopies.filtered-intersect | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $d/expressed_retrocopies.bed

done
