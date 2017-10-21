#!/bin/bash

INPUT=$1
OUTPUT=$2

cat $INPUT | awk '{print $1,$2,$3,$5,0,$4}' | sed -r 's/\s+/\t/g' | sort -k1,1 -k2,2n > $OUTPUT

