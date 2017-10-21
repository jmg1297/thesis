#!/bin/bash

SRCDIR=data/all_exons_int_rt_blocks
TARGETDIR=data/rt_transcripts_vs_ensembl
ENSEMBLBED=$HOME/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.protein_coding_transcripts.sorted.bed

for d in $(find $SRCDIR -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  echo $MRG
  mkdir $TARGETDIR/$MRG
  bedtools intersect -wao -a $d/transcripts_rt_overlap_summ.bed -b $ENSEMBLBED > $TARGETDIR/$MRG/rt_transcripts-ensembl_prot_coding.intersect
  bedtools intersect -wao -a $d/transcripts_rt_overlap_summ.filter_50pc.bed -b $ENSEMBLBED > $TARGETDIR/$MRG/rt_transcripts_gt50pc-ensembl_prot_coding.intersect
done
