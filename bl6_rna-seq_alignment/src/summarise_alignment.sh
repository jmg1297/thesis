#!/bin/bash

python src/summarise_alignment.py data/alignment_summary.json data/star_aligned data/rrna_filtered alignment_summ.png
mv *alignment_summ.png data/alignment_summaries/
