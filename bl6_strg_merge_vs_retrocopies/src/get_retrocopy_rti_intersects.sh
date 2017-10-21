#!/bin/bash

RTIBED=$HOME/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed
ALLRCBED=$HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed
TARGETDIR=data/ALL_expr_retrocopies_vs_RTIs

bedtools intersect -wo -f 0.9 -a data/all_exons_int_retrocopy_blocks/ALL/expressed_retrocopies.bed -b $RTIBED > $TARGETDIR/ALL_expr_rc-rtis.intersect

bedtools intersect -wo -f 0.9 -a $ALLRCBED -b $RTIBED > $TARGETDIR/all_rc-rtis.intersect

cat $TARGETDIR/ALL_expr_rc-rtis.intersect | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $TARGETDIR/ALL_expr_rcs_in_rtis.bed

cat $TARGETDIR/all_rc-rtis.intersect | awk '{print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g' > $TARGETDIR/all_rcs_in_rtis.bed


