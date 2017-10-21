#!/bin/bash

MBEDFNAME=data/exons_int_rt_blocks/transcripts_retrotransposons.mbed
TARGETFNAME=data/exons_int_rt_blocks/rts_intersecting_transcripts

cat $MBEDFNAME | awk '{print $7,$8,$9,$10,$11,$12}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $TARGETFNAME".bed"
cat $TARGETFNAME".bed" | grep "LINE" > $TARGETFNAME."LINE.bed"
cat $TARGETFNAME".bed" | grep "SINE" > $TARGETFNAME."SINE.bed"
cat $TARGETFNAME".bed" | grep "LTR" > $TARGETFNAME."LTR.bed"
