#!/bin/bash

GTF=data/stringtie_merge/ENCFF001RT_merge.gtf

cat $GTF | \
grep -P "\sexon\s" | \
awk '{print $1,$4,$5,$12,$14,0,$7}' | \
sed -r -e 's/[\";]//g' -e 's/ /./4' -e 's/\s+/\t/g' -e 's/^/chr/' | \
sort -k1,1 -k2,2n \
> data/exons/exons.bed
