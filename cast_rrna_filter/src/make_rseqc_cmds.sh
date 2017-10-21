#!/bin/bash

TARGETDIR=data/rrna_filter
RRNA_BED=$HOME/projects/cast_rrna_identification/data/blast_rRNA_to_cast/cast_rRNA.bed
SRCDIR=$HOME/projects/cast_rna-seq_transcriptome_reconstruction/data/star_aligned

for d in $(find $SRCDIR -type d -name "*CAST*");
do
  SAMPLE=$(echo $d | cut -d/ -f8)
  echo "mkdir "$TARGETDIR"/"$SAMPLE"; python "$HOME"/bin/split_bam.py -i "$d"/Aligned.out.bam -r "$RRNA_BED" -o "$TARGETDIR"/"$SAMPLE"/rseqc_rrna" >> src/rseqc_rrna.cmds
done
