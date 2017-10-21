#!/bin/bash

/usr/bin/python /home/jg600/projects/stringtie_merge/src/make_strg_merge_cmds.py \
/home/jg600/projects/data/rna-seq.fastq/metadata.db \
CAST-Eij \
0 \
/home/jg600/projects/cast_stringtie/data/stringtie/ \
out.gtf \
data/stringtie_merge/ \
src/stringtie_merge.cmds

