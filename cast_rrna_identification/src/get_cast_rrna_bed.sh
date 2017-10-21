#!/bin/bash

cd data/blast_rRNA_to_cast

bedtools merge \
-i <(cat mm10_rRNA-CAST_genome.blast.out.plusonly | awk '{print $2,$9,$10,$1}' | sort -k1,1 -k2,2n | uniq | sed -r 's/\s+/\t/g') \
-c 4 -o distinct | \
awk '{print $1,$2,$3,$4,0,".",$2,$3,0,1,$3-$2,0}' | sed -r 's/\s+/\t/g' \
> cast_rRNA.bed
