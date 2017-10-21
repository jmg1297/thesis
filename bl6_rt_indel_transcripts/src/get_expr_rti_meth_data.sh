#!/bin/bash

SCRIPTDIR=$HOME/projects/strg_merge_all_exons_vs_retrocopies/src
BEDGRAPH=/home/jg600/bs-seq.bedgraphs/merged/BL6_complete.bedgraph

python $SCRIPTDIR/get_meth_data.py \
-i data/expr_rti_contents/ALL/expressed_rtis.classified \
-l EXPRESSED_RTIs \
-i $HOME/projects/mm10_ref/repeatmasker/retrotransposon_indels/rt_indels.numbered.bed \
-l ALL_RTIs \
-g $BEDGRAPH \
-L BL6MRG \
6 data/methylation/ALL_expr_rtis_vs_all_rtis_methylation_data.json


