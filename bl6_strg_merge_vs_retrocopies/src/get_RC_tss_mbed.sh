#!/bin/bash

MBED=data/all_exons_int_retrocopy_blocks/ALL/transcript_retrocopy_parent.mbed
TARGET=data/ALL_rc_transcript_tss/tss_retrocopy_parent.mbed

cat $MBED | perl -lane 'if ($F[5] eq "-"){$F[1]=$F[2];} else {$F[2]=$F[1]}; print "@F"' | sed -r 's/\s+/\t/g' > $TARGET
