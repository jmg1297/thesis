#!/bin/bash

SCRIPTDIR=$HOME/thesis_projects/src
BTDIR=$HOME/thesis_projects/strg_merge_all_exons_vs_rt_blocks/data/all_exons_int_rt_blocks
LIVERDIR=data/exons_int_rt_blocks

python $SCRIPTDIR/venn_diagrams.py \
-i $LIVERDIR/rts_intersecting_transcripts.bed -l Liver \
-i $BTDIR/B/rts_intersecting_transcripts.bed -l B \
-i $BTDIR/T/rts_intersecting_transcripts.bed -l T \
imgs/venn_BTLiver_rts_intersecting_transcripts.svg

python $SCRIPTDIR/venn_diagrams.py \
-i $LIVERDIR/rts_intersecting_gt50pc_transcripts.bed -l Liver \
-i $BTDIR/B/rts_intersecting_gt50pc_transcripts.bed -l B \
-i $BTDIR/T/rts_intersecting_gt50pc_transcripts.bed -l T \
imgs/venn_BTLiver_rts_intersecting_gt50pc_transcripts.svg

python $SCRIPTDIR/venn_diagrams.py \
-i $LIVERDIR/rts_intersecting_gt90pc_transcripts.bed -l Liver \
-i $BTDIR/B/rts_intersecting_gt90pc_transcripts.bed -l B \
-i $BTDIR/T/rts_intersecting_gt90pc_transcripts.bed -l T \
imgs/venn_BTLiver_rts_intersecting_gt90pc_transcripts.svg
