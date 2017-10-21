#!/bin/bash

MBED=data/all_exons_int_retrocopy_blocks/ALL/transcript_retrocopy_parent.mbed
TSSBED=data/ALL_rc_transcript_tss/rc_transcript_tss.bed

cat $MBED | \
awk '{print $1,$2,$3,$4,$5,$6}' | \
perl -lane 'if ($F[5] eq "+"){$F[2]=$F[1]} else {$F[1]=$F[2]}; print "@F"' | \
sed -r 's/\s+/\t/g' \
> $TSSBED
