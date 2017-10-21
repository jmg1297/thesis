#!/bin/bash

INPUT=$1
OUTPUT=$2

cat $INPUT |\
awk '{print $2,$9,$10,$1,0}' | \
perl -lane 'my $strand; my $tmp; if ($F[2] >= $F[1]){$strand="+"} else {$strand="-"; $tmp=$F[2]; $F[2]=$F[1]; $F[1]=$tmp}; print "@F\t$strand"' |\
sort -k1,1 -k2,2n |\
sed -r 's/\s+/\t/g' \
> $OUTPUT
