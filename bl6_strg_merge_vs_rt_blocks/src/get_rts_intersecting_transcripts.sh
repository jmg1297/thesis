#!/bin/bash

SRCDIR=data/all_exons_int_rt_blocks
MBEDFNAME=transcripts_retrotransposons.mbed
TARGETFNAME=rts_intersecting_transcripts

for d in $(find $SRCDIR -type d | sed '1d');
do
  cat $d/$MBEDFNAME | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $d/$TARGETFNAME".bed"
  cat $d/$TARGETFNAME".bed" | grep "LINE" > $d/$TARGETFNAME."LINE.bed"
  cat $d/$TARGETFNAME".bed" | grep "SINE" > $d/$TARGETFNAME."SINE.bed"
  cat $d/$TARGETFNAME".bed" | grep "LTR" > $d/$TARGETFNAME."LTR.bed"
done
