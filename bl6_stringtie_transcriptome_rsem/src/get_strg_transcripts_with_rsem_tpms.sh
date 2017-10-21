#!/bin/bash

STRGDIR=$HOME/projects/rna-seq_transcriptome_reconstruction/data/stringtie
RSEMDIR=$HOME/projects/stringtie_transcriptome_rsem/data/stringtie_rsem_expression
OUTDIR=$HOME/projects/stringtie_transcriptome_rsem/data/strg_transcripts_rsem_tpm

for f in $(find $STRGDIR -type d -name "*C57BL-6J*");
do
SAMPLE=$(echo $f | cut -d/ -f8)
STRG=$STRGDIR/$SAMPLE/transcripts.bed
RSEM=$RSEMDIR/$SAMPLE/stringtie.isoforms.results
OUT=$OUTDIR/$SAMPLE/transcripts.bed
echo $SAMPLE
mkdir $OUTDIR/$SAMPLE; paste <(cat $STRG | grep -vP "^[MGJ]" | awk '{print $1,$2,$3,$4,$6}' | LC_ALL=C sort -k4) <(cat $RSEM | awk '{print $1,$6}' | sed '1d' | LC_ALL=C sort -k1) | awk '{print $1,$2,$3,$4,$7,$5}' | sed -r 's/\s+/\t/g' | sort -k1,1 -k2,2n > $OUT
done

