#!/bin/bash

for d in $(find data/split_by_gffcmp_class -type d | sed '1d');
do
MRG=$(echo $d | cut -d/ -f3)
mkdir data/novel_exons/$MRG
cat $d/novel_non-coding.gtf | grep -P "\sexon\s" | awk '{print $1,$4,$5,$12,$14,0,$7}' | sed -r -e 's/[\";]//g' -e 's/ /./4' -e 's/^/chr/g' -e 's/\s+/\t/g' | sort -k1,1 -k2,2n > data/novel_exons/$MRG/novel_exons.bed
done
