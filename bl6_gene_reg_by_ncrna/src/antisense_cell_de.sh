#!/bin/bash

# First get the list of transcripts in the introns of a reference transcript
# for B cells and T cells

TMAPDIR=$HOME/thesis_projects/stringtie_merge/data/gffcompare_ensembl
TARGETDIR=data/with_antisense
BALLGOWN=$HOME/thesis_projects/ensembl_ballgown/data/ballgown/cell_de_results_transcripts.tsv

cat $TMAPDIR/B/B.B_merge.gtf.tmap | grep -P "\sx\s" | awk '{print $2}' | sort | uniq > $TARGETDIR/B_ensembl_with_antisense.txt
cat $TMAPDIR/T/T.T_merge.gtf.tmap | grep -P "\sx\s" | awk '{print $2}' | sort | uniq > $TARGETDIR/T_ensembl_with_antisense.txt

comm -12 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/cell_shared.ensembl_with_antisense.txt
comm -13 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/T_spec.ensembl_with_antisense.txt
comm -23 $TARGETDIR/B_ensembl_with_antisense.txt $TARGETDIR/T_ensembl_with_antisense.txt > $TARGETDIR/B_spec.ensembl_with_antisense.txt

python src/ensembl_fc_volcano.py \
$TARGETDIR/B_spec.ensembl_with_antisense.txt \
$TARGETDIR/T_spec.ensembl_with_antisense.txt \
$TARGETDIR/cell_shared.ensembl_with_antisense.txt \
$BALLGOWN \
imgs/volcano_ensembl_with_antisense_cell_de.svg ""
