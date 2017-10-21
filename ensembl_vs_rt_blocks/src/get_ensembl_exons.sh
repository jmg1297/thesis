#!/bin/bash

GTF=/home/jg600/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf
EXONS=/home/jg600/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.exons.bed

cat $GTF | \
sed '/^#/d' | \
grep -P "\sexon\s" | \
awk '{print $1,$4,$5,$14,$18,"0",$7,$24}' | \
sed -r -e 's/"; "/./' -e 's/[";]//g' -e 's/^/chr/' -e 's/\s+/\t/g' | \
sort -k1,1 -k2,2n \
> $EXONS
