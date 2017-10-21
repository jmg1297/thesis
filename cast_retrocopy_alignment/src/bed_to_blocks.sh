#!/bin/bash

INPUT=$1
OUTPUT=$2

cat $INPUT | \
perl -lane '$F[3] = "retrocopy:putative:$F[3]:$..1"; print "@F"' | \
sed -r 's/\s+/\t/g' \
> $OUTPUT
