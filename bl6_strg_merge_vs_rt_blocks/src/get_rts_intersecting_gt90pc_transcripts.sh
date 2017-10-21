#!/bin/bash

SRCDIR=data/all_exons_int_rt_blocks
SUMMFNAME=transcripts_rt_overlap_summ.filter_90pc.bed
MBEDFNAME=transcripts_retrotransposons.mbed
TARGETFNAME=rts_intersecting_gt90pc_transcripts

for d in $(find $SRCDIR -type d | sed '1d');
do
  grep \
    -f <(cat $d/$SUMMFNAME | awk '{print $1,$2,$3,$4,$5,$6}' | sed -r 's/\s+/\t/g') \
    $d/$MBEDFNAME | \
  awk '{print $7,$8,$9,$10,$11,$12}' | \
  sort -k1,1 -k2,2n | \
  uniq | \
  sed -r 's/\s+/\t/g' \
  > $d/$TARGETFNAME".bed"

  cat $d/$TARGETFNAME".bed" | grep "LINE" > $d/$TARGETFNAME".LINE.bed"
  cat $d/$TARGETFNAME".bed" | grep "SINE" > $d/$TARGETFNAME".SINE.bed"
  cat $d/$TARGETFNAME".bed" | grep "LTR" > $d/$TARGETFNAME".LTR.bed"
done
