#!/bin/bash

CHIPSEQ=$HOME/projects/data/bp_histone_chip-seq
H3K27ac=$CHIPSEQ/H3K27ac/C57BL-6J/merged/C57BL-6J_H3K27ac.bed
H3K27me3=$CHIPSEQ/H3K27me3/C57BL-6J/merged/C57BL-6J_H3K27me3.bed
H3K36me3=$CHIPSEQ/H3K36me3/C57BL-6J/merged/C57BL-6J_H3K36me3.bed
H3K4me3=$CHIPSEQ/H3K4me3/C57BL-6J/merged/C57BL-6J_Male_H3K4me3.bed
ALLRTIBED=$HOME/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed
SRCDIR=data/expr_rti_contents/ALL

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.classified \
$ALLRTIBED \
500 500 16 \
"Histone Modification Peaks Near Expressed RTIs vs All RTIs" \
imgs/boxplot_hmods_in_rti_windows.svg

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.LINE_classified \
$ALLRTIBED \
35 500 16 \
"Histone Modification Peaks Near Expressed LINE RTIs vs All RTIs" \
imgs/boxplot_hmods_in_LINE_rti_windows.svg

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.SINE_classified \
$ALLRTIBED \
65 500 16 \
"Histone Modification Peaks Near Expressed SINE RTIs vs All RTIs" \
imgs/boxplot_hmods_in_SINE_rti_windows.svg

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.LTR_classified \
$ALLRTIBED \
85 500 16 \
"Histone Modification Peaks Near Expressed LTR RTIs vs All RTIs" \
imgs/boxplot_hmods_in_LTR_rti_windows.svg

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.pseudogene_classified \
$ALLRTIBED \
255 500 16 \
"Histone Modification Peaks Near Expressed Retrocopy RTIs vs All RTIs" \
imgs/boxplot_hmods_in_retrocopy_rti_windows.svg

python src/boxplots_histone_mod_windows_from_rtis.py \
-b $H3K27ac -l H3K27ac \
-b $H3K27me3 -l H3K27me3 \
-b $H3K36me3 -l H3K36me3 \
-b $H3K4me3 -l H3K4me3 \
-d 0 -d 500 -d 1000 -d 1500 -d 2000 -d 2500 -d 5000 \
$SRCDIR/expressed_rtis.NONE_classified \
$ALLRTIBED \
115 500 16 \
"Histone Modification Peaks Near Expressed Unclassified RTIs vs All RTIs" \
imgs/boxplot_hmods_in_NONE_rti_windows.svg
