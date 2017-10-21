#!/bin/bash

HMODDIR=$HOME/projects/data/bp_histone_chip-seq

python src/plot_transcript_hmods.py \
-b $HMODDIR/H3K27ac/C57BL-6J/merged/C57BL-6J_B_H3K27ac.bed -g B -m H3K27ac \
-b $HMODDIR/H3K27ac/C57BL-6J/merged/C57BL-6J_T_H3K27ac.bed -g T -m H3K27ac \
-b $HMODDIR/H3K27me3/C57BL-6J/merged/C57BL-6J_B_H3K27me3.bed -g B -m H3K27me3 \
-b $HMODDIR/H3K27me3/C57BL-6J/merged/C57BL-6J_T_H3K27me3.bed -g T -m H3K27me3 \
-b $HMODDIR/H3K36me3/C57BL-6J/merged/C57BL-6J_Male_B_H3K36me3.bed -g B -m H3K36me3 \
-b $HMODDIR/H3K36me3/C57BL-6J/merged/C57BL-6J_Male_T_H3K36me3.bed -g T -m H3K36me3 \
-b $HMODDIR/H3K4me3/C57BL-6J/merged/C57BL-6J_Male_B_H3K4me3.bed -g B -m H3K4me3 \
-b $HMODDIR/H3K4me3/C57BL-6J/merged/C57BL-6J_Male_T_H3K4me3.bed -g T -m H3K4me3 \
-l 5000 -r 5000 \
data/all_exons_int_retrocopy_blocks/T_spec_parents_with_T_upreg.txt \
/home/jg600/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf \
8 imgs/T_spec_parents_T_upreg_hmods/
