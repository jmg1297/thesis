#!/bin/bash

SRCDIR=data/exons_int_rt_blocks
SUMMFNAME=transcripts_rt_overlap_summ.filter_90pc.bed
MBEDFNAME=transcripts_retrotransposons.mbed
TARGETFNAME=$SRCDIR/rts_intersecting_gt90pc_transcripts

grep \
  -f <(cat $SRCDIR/$SUMMFNAME | awk '{print $1,$2,$3,$4,$5,$6}' | sed -r 's/\s+/\t/g') \
  $SRCDIR/$MBEDFNAME | \
awk '{print $7,$8,$9,$10,$11,$12}' | \
sort -k1,1 -k2,2n | \
uniq | \
sed -r 's/\s+/\t/g' \
> $TARGETFNAME".bed"

cat $TARGETFNAME".bed" | grep "LINE" > $TARGETFNAME".LINE.bed"
cat $TARGETFNAME".bed" | grep "SINE" > $TARGETFNAME".SINE.bed"
cat $TARGETFNAME".bed" | grep "LTR" > $TARGETFNAME".LTR.bed"
