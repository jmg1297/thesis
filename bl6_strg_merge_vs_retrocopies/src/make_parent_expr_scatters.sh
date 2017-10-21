#!/bin/bash

/usr/bin/python src/scatter_cell_spec_parent_expr_by_rel_sense.py \
data/all_exons_int_retrocopy_blocks/B/spec_expr_retrocopies.mbed \
data/all_exons_int_retrocopy_blocks/T/spec_expr_retrocopies.mbed \
/home/jg600/projects/kallisto_ensembl/data/kallisto_results.db \
data/strand_combos.json \
0 \
1 \
$'Parent Expression for \nCell-Type-Specifically-Expressed Retrocopies\n' \
imgs/cell_spec_parent_expr_scatter_by_rel_sense.svg

/usr/bin/python src/scatter_sex_spec_parent_expr_by_rel_sense.py \
data/all_exons_int_retrocopy_blocks/Female/spec_expr_retrocopies.mbed \
data/all_exons_int_retrocopy_blocks/Male/spec_expr_retrocopies.mbed \
/home/jg600/projects/kallisto_ensembl/data/kallisto_results.db \
data/strand_combos.json \
0 \
1 \
$'Parent Expression for \nSex-Specifically-Expressed Retrocopies' \
imgs/sex_spec_parent_expr_scatter_by_rel_sense.svg
