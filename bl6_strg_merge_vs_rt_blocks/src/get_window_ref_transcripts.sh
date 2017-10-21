#!/bin/bash

WINDOW_SIZES=(1000 5000 10000 50000)
TRANSCRIPTS=$HOME/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.transcripts.bed
SRCDIR=data/all_exons_int_rt_blocks
TARGETDIR=data/ref_transcript_expr_near_rts

for W in ${WINDOW_SIZES[@]};
do
  echo $W
  bedtools window -w $W -a $SRCDIR/B/rts_intersecting_transcripts.bed -b $TRANSCRIPTS > $TARGETDIR/B_rts-ensembl.$W.window
  bedtools window -w $W -a $SRCDIR/T/rts_intersecting_transcripts.bed -b $TRANSCRIPTS > $TARGETDIR/T_rts-ensembl.$W.window

  cat $TARGETDIR/B_rts-ensembl.$W.window | awk '{print $10}' | sort | uniq > $TARGETDIR/B_rts-ensembl_ids-$W-window.txt
  cat $TARGETDIR/T_rts-ensembl.$W.window | awk '{print $10}' | sort | uniq > $TARGETDIR/T_rts-ensembl_ids-$W-window.txt

  comm -12 $TARGETDIR/B_rts-ensembl_ids-$W-window.txt $TARGETDIR/T_rts-ensembl_ids-$W-window.txt > $TARGETDIR/cell_shared.rts-ensembl_ids-$W-window.txt
  comm -23 $TARGETDIR/B_rts-ensembl_ids-$W-window.txt $TARGETDIR/T_rts-ensembl_ids-$W-window.txt > $TARGETDIR/B_spec.rts-ensembl_ids-$W-window.txt
  comm -13 $TARGETDIR/B_rts-ensembl_ids-$W-window.txt $TARGETDIR/T_rts-ensembl_ids-$W-window.txt > $TARGETDIR/T_spec.rts-ensembl_ids-$W-window.txt
done

bedtools intersect -wo -a $SRCDIR/B/rts_intersecting_transcripts.bed -b $TRANSCRIPTS > $TARGETDIR/B_rts-ensembl.intersect
bedtools intersect -wo -a $SRCDIR/T/rts_intersecting_transcripts.bed -b $TRANSCRIPTS > $TARGETDIR/T_rts-ensembl.intersect

cat $TARGETDIR/B_rts-ensembl.intersect | awk '{print $10}' | sort | uniq > $TARGETDIR/B_rts-ensembl_ids-intersect.txt
cat $TARGETDIR/T_rts-ensembl.intersect | awk '{print $10}' | sort | uniq > $TARGETDIR/T_rts-ensembl_ids-intersect.txt

comm -12 $TARGETDIR/B_rts-ensembl_ids-intersect.txt $TARGETDIR/T_rts-ensembl_ids-intersect.txt > $TARGETDIR/cell_shared.rts-ensembl_ids-intersect.txt
comm -23 $TARGETDIR/B_rts-ensembl_ids-intersect.txt $TARGETDIR/T_rts-ensembl_ids-intersect.txt > $TARGETDIR/B_spec.rts-ensembl_ids-intersect.txt
comm -13 $TARGETDIR/B_rts-ensembl_ids-intersect.txt $TARGETDIR/T_rts-ensembl_ids-intersect.txt > $TARGETDIR/T_spec.rts-ensembl_ids-intersect.txt
