#!/bin/bash

MBED=data/all_exons_int_retrocopy_blocks/ALL/transcript_retrocopy_parent.mbed

cat $MBED | perl -lane 'BEGIN{my $tss;} if ($F[5] eq "+"){$tss=$F[1]} else {$tss=$F[2]}; if ($F[7] < $tss & $tss < $F[8] & $F[5] eq $F[11]){print $_}' > data/ALL_rc_transcript_tss/tss_inside_rc.mbed

