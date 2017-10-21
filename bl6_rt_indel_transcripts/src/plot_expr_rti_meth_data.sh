#!/bin/bash

SCRIPTDIR=$HOME/projects/strg_merge_all_exons_vs_retrocopies/src

python $SCRIPTDIR/plot_meth_data.py data/methylation/ALL_expr_rtis_vs_all_rtis_methylation_data.json "Methylation of expressed RTIs vs All RTIs" imgs/expr_rti_vs_all_rti_meth.svg
