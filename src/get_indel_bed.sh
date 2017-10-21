#!/bin/bash

zcat cast_ref/genome/downloads/CAST_EiJ.mgp.v5.indels.dbSNP142.normed.vcf.gz | \
sed '/^#/d' | \
awk '{print $1,$2,$4,$5}' | \
perl -lane 'BEGIN{my $len; my $end;} $len = length $F[2]; $end = $F[1]+$len-1; print join("\t", ($F[0], $F[1], $end, "@F[2..$#F]"))' | \
awk '{print $1,$2,$3,"CL=INDEL;REF="$4";ALT="$5}' | \
sed -r -e 's/^/chr/' -e 's/\s+/\t/g' | \
sort -k1,1 -k2,2n \
> cast_ref/genome/variation/CAST_EiJ.mgp.v5.indels.dbSNP142.bed
