#!/bin/bash

python /home/jg600/projects/rna-seq_transcriptome_reconstruction/src/summarise_alignment.py \
data/alignment_summary.json \
data/star_aligned/ \
/home/jg600/projects/cast_rrna_filter/data/rrna_filter/ \
alignment_summary.svg

mv *alignment_summary.svg imgs
