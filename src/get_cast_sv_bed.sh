#!/bin/bash

zcat cast_ref/genome/downloads/28strains.REL-1410-SV.sdp.tab.gz | \
sed '1,6d' | \
awk '{print $4,$13}' | \
sed -r -e '/\s+0$/d' -e 's/;/ /' -e '1d' | \
awk '{print $2,$1,$3}' | \
sed -r 's/[:;-]/ /g' | \
awk '{print $1,$2,$3,$5,$8,$6,$9,$7,$10}' | \
sed -r -e 's/ /=/8' -e 's/ /=/6' -e 's/ /=/4' -e 's/ /,/4g' -e 's/\s+/\t/g' | \
sort -k1,1 -k2,2n | \
sed -r 's/^/chr/' \
> cast_ref/genome/variation/CAST-Eij.REL-1410-SV.sdp.bed
