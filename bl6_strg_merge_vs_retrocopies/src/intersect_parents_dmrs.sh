#!/bin/bash

RCDIR=data/all_exons_int_retrocopy_blocks
LMRDIR=/home/nw11/data/blueprint/wgbs-t-b-cell-bl6-lmrs
TARGETDIR=data/parent_methylation
CHROMSIZES=$HOME/projects/mm10_ref/genome/mm10.chrom.sizes

# First, get a BED file of parents with added sequence either side
cat $RCDIR/B/spec_expr_retrocopies.mbed | cut -f14-19 | sort -k1,1 -k2,2n | uniq > $TARGETDIR/B_spec_parents.bed
cat $RCDIR/T/spec_expr_retrocopies.mbed | cut -f14-19 | sort -k1,1 -k2,2n | uniq > $TARGETDIR/T_spec_parents.bed

bedtools slop -i $TARGETDIR/B_spec_parents.bed -h $CHROMSIZES -b 1000 > $TARGETDIR/B_spec_parents.slop1000b.bed
bedtools slop -i $TARGETDIR/T_spec_parents.bed -h $CHROMSIZES -b 1000 > $TARGETDIR/T_spec_parents.slop1000b.bed

# Now get BED files for the LMRs

