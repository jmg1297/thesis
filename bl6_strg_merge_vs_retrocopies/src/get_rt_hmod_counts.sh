#!/bin/bash

HMODBEDS=( \
$HOME/projects/data/bp_histone_chip-seq/H3K27ac/C57BL-6J/merged/C57BL-6J_B_H3K27ac.bed \
$HOME/projects/data/bp_histone_chip-seq/H3K27me3/C57BL-6J/merged/C57BL-6J_B_H3K27me3.bed \
$HOME/projects/data/bp_histone_chip-seq/H3K36me3/C57BL-6J/merged/C57BL-6J_H3K36me3.bed \
$HOME/projects/data/bp_histone_chip-seq/H3K4me3/C57BL-6J/merged/C57BL-6J_Male_H3K4me3.bed
)

# Only do ALL merge, for now
# First do the unexpressed retrocopies with randomly assigned strands
BASEDIR=$(pwd)
OUTFILE=$BASEDIR/data/expr_retrocopies_vs_retrotransposon_blocks/ALL/hmod_overlap_counts.txt

python $HOME/projects/src/get_random_regions.py \
$HOME/projects/mm10_ref/genome/mm10.chrom.sizes \
$HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed \
50 \
18000 \
data/tmp.bed

sort -k1,1 -k2,2n data/tmp.bed > data/expr_retrocopies_vs_retrotransposon_blocks/ALL/rand_rclike_regions.18000.bed
rm data/tmp.bed

echo "HERE"

bedtools closest -D a -id -io \
-a data/expr_retrocopies_vs_retrotransposon_blocks/ALL/rand_rclike_regions.18000.bed \
-b $HOME/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.blocks.formatted.bed \
> data/expr_retrocopies_vs_retrotransposon_blocks/ALL/rand_rclike_regions.18000.bed.rt_blocks.closest

echo "NOW HERE"

RANDCLOSEST=$BASEDIR/data/expr_retrocopies_vs_retrotransposon_blocks/ALL/rand_rclike_regions.18000.bed.rt_blocks.closest

cd data/unexpr_retrocopies_with_rand_strand/ALL

for h in ${HMODBEDS[@]};
do
  cd $BASEDIR/data/unexpr_retrocopies_with_rand_strand/ALL
  for f in $(ls *.closest);
  do
    COUNT=$(bedtools intersect -u \
    -a <(cat $f | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq) \
    -b $h | \
    wc -l)

    TOTAL=$(cat $f | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq | wc -l)

    HMOD=$(echo $h | cut -d/ -f7)

    echo $HMOD","$f","$COUNT","$TOTAL >> $OUTFILE
  done

  cd $BASEDIR/data/expr_retrocopies_vs_retrotransposon_blocks/ALL
  F=retrocopy_transcripts-rt_blocks.closest
  COUNT=$(bedtools intersect -u \
  -a <(cat $F | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq) \
  -b $h | \
  wc -l)
  TOTAL=$(cat $F | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq | wc -l)
  echo $HMOD","$F","$COUNT","$TOTAL >> $OUTFILE

  COUNT=$(bedtools intersect -u \
  -a <(cat $RANDCLOSEST | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq) \
  -b $h | \
  wc -l)
  TOTAL=$(cat $RANDCLOSEST | cut -f7-12 | sed '/\s\.\s/d' | sort -k1,1 -k2,2n | uniq | wc -l)
  echo $HMOD","$RANDCLOSEST","$COUNT","$TOTAL >> $OUTFILE

done
