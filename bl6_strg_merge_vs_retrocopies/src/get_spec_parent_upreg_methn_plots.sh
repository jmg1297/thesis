#!/bin/bash

METHDIR=$HOME/bs-seq.bedgraphs/smoothed

python src/plot_transcript_methylation.py \
-b $METHDIR/BL6-F-B-108.smoothed.bedgraph -g B -L BL6-F-B-108 \
-b $METHDIR/BL6-F-B-138.smoothed.bedgraph -g B -L BL6-F-B-138 \
-b $METHDIR/BL6-M-B-88.smoothed.bedgraph -g B -L BL6-M-B-88 \
-b $METHDIR/BL6-M-B-98.smoothed.bedgraph -g B -L BL6-M-B-98 \
-b $METHDIR/BL6-F-T-7905947.smoothed.bedgraph -g T -L BL6-F-T-7905947 \
-b $METHDIR/BL6-F-T-7905948.smoothed.bedgraph -g T -L BL6-F-T-7905948 \
-b $METHDIR/BL6-M-T-7905953.smoothed.bedgraph -g T -L BL6-M-T-7905953 \
-b $METHDIR/BL6-M-T-7905954.smoothed.bedgraph -g T -L BL6-M-T-7905954 \
data/all_exons_int_retrocopy_blocks/B_spec_parents_with_B_upreg.txt \
/home/jg600/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf \
4 imgs/B_spec_parents_B_upreg_methn


python src/plot_transcript_methylation.py \
-b $METHDIR/BL6-F-B-108.smoothed.bedgraph -g B -L BL6-F-B-108 \
-b $METHDIR/BL6-F-B-138.smoothed.bedgraph -g B -L BL6-F-B-138 \
-b $METHDIR/BL6-M-B-88.smoothed.bedgraph -g B -L BL6-M-B-88 \
-b $METHDIR/BL6-M-B-98.smoothed.bedgraph -g B -L BL6-M-B-98 \
-b $METHDIR/BL6-F-T-7905947.smoothed.bedgraph -g T -L BL6-F-T-7905947 \
-b $METHDIR/BL6-F-T-7905948.smoothed.bedgraph -g T -L BL6-F-T-7905948 \
-b $METHDIR/BL6-M-T-7905953.smoothed.bedgraph -g T -L BL6-M-T-7905953 \
-b $METHDIR/BL6-M-T-7905954.smoothed.bedgraph -g T -L BL6-M-T-7905954 \
data/all_exons_int_retrocopy_blocks/T_spec_parents_with_T_upreg.txt \
/home/jg600/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf \
4 imgs/T_spec_parents_T_upreg_methn
