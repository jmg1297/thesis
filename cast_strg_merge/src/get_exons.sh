#!/bin/bash

for f in $(find data/stringtie_merge -type f -name "*_merge.gtf");
do
  MRG=$(echo $f | cut -d/ -f3 | sed -r 's/_merge\.gtf//')
  mkdir data/exons/$MRG
  cat $f | grep -P "\sexon\s" | awk '{print $1,$4,$5,$12,$14,0,$7}' | perl -lane 'if ($F[0] !~ /^chr/){;} else {print $_}' | sed -r -e 's/[\";]//g' -e 's/ /./4' -e 's/\s+/\t/g' | sort -k1,1 -k2,2n > data/exons/$MRG/exons.bed
done
