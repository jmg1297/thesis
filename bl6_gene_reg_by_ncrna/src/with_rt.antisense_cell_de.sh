#!/bin/bash

# First get the list of transcripts in the introns of a reference transcript
# for B cells and T cells

BEDDIR=$HOME/thesis_projects/strg_merge_all_exons_vs_rt_blocks/data/all_exons_int_rt_blocks
TARGETDIR=data/with_rt_antisense
BALLGOWN=$HOME/thesis_projects/ensembl_ballgown/data/ballgown/cell_de_results_transcripts.tsv

cat $BEDDIR/B/transcripts_rt_overlap_summ+gffcmp.bed | grep -P "\sx\s" | awk '{print $9}' | sort | uniq > $TARGETDIR/B_ensembl_with_antisense.txt
cat $BEDDIR/T/transcripts_rt_overlap_summ+gffcmp.bed | grep -P "\sx\s" | awk '{print $9}' | sort | uniq > $TARGETDIR/T_ensembl_with_antisense.txt

comm -12 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/cell_shared.ensembl_with_antisense.txt
comm -13 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/T_spec.ensembl_with_antisense.txt
comm -23 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/B_spec.ensembl_with_antisense.txt

python src/ensembl_fc_volcano.py \
$TARGETDIR/B_spec.ensembl_with_antisense.txt \
$TARGETDIR/T_spec.ensembl_with_antisense.txt \
$TARGETDIR/cell_shared.ensembl_with_antisense.txt \
$BALLGOWN \
imgs/volcano_ensembl_with_rt_antisense_cell_de.svg ""
