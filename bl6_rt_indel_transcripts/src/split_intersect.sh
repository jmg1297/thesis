#!/bin/bash

cd data/strg_merge_int_rt_indels

for d in $(find . -type d | sed '1d');
do
  cd $d
  cat transcripts.rt_indels.intersect | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > rti_transcripts.bed
  cat transcripts.rt_indels.intersect | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > expressed_rtis.bed
  cd ..
done
