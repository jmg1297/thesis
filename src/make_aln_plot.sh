#!/bin/bash

/usr/bin/python src/density_plot_align_scores.py \
-i strg_merge_exons_vs_retrocopies/data/novel_exons_int_retrocopy_blocks/ALL/expressed_retrocopies_sense_rel_parent.bed -l SENSE_EXPR -c blue \
-i strg_merge_exons_vs_retrocopies/data/novel_exons_int_retrocopy_blocks/ALL/expressed_retrocopies_antisense_rel_parent.bed -l ANTISENSE_EXPR -c red \
-i strg_merge_exons_vs_retrocopies/data/novel_exons_int_retrocopy_blocks/ALL/expressed_retrocopies.bed -l ALL_EXPR -c orange \
-i strg_merge_exons_vs_retrocopies/data/novel_exons_int_retrocopy_blocks/ALL/unexpressed_retrocopies.bed -l UNEXPR -c green \
-i mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed  -l ALL_RETROCOPIES -c purple \
mm10_ref/retrogenes/parent_alignments/ \
mm10_ref/retrogenes/random_alignments/ \
mm10_ref/retrogenes/retrogene_parent.db \
"Retrocopy/Parent Alignment Scores" \
strg_merge_exons_vs_retrocopies/imgs/aln_score_density_plot_ALL.svg
