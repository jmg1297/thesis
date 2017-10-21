#!/bin/bash

# First use python to get sets of samples from
# the metadata database and write the commands
# to a text file, then run the commands using
# parallel.

# Use the system python install for sqlite3 module
/usr/bin/python src/make_strg_merge_cmds.py \
/home/jg600/projects/data/rna-seq.fastq/metadata.db  \
/home/jg600/projects/rna-seq_transcriptome_reconstruction/data/stringtie \
out.gtf \
data/stringtie_merge \
src/strg_merge.cmds

cat src/strg_merge.cmds | parallel -j $1
