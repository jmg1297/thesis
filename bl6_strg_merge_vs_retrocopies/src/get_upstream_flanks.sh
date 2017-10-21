#!/bin/bash

FLANKSIZE=$1

BASEDIR=$(pwd)
TARGETDIR=data/expr_retrocopies_upstream_flanks
GENOME=$HOME/projects/mm10_ref/genome/mm10.chrom.sizes
UNEXPRDIR=data/unexpr_retrocopies_with_rand_strand

# First do the retrocopy transcripts

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  mkdir $TARGETDIR/$MRG
  bedtools flank -i $d/transcripts.retrocopy.bed -g $GENOME -s -l $FLANKSIZE -r 0 > $TARGETDIR/$MRG/retrocopy_transcripts.flank_sl1000r0.bed

  # Now do the unexpressed retrocopies with random strands
  cd $UNEXPRDIR/$MRG
  for f in $(find . -type f -name "*.bed");
  do
    bedtools flank -i $f -g $GENOME -s -l $FLANKSIZE -r 0 > $f".flank_sl1000r0.bed"
  done
  cd $BASEDIR
done

# Finally generate some retrocopy-like regions and get their flanks

RANDBED=data/expr_retrocopies_upstream_flanks/ALL/rand_rclike_regions.18000.bed

python $HOME/projects/src/get_random_regions.py \
$HOME/projects/mm10_ref/genome/mm10.chrom.sizes \
$HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed \
50 \
18000 \
data/tmp.bed

sort -k1,1 -k2,2n data/tmp.bed > $RANDBED
rm data/tmp.bed

bedtools flank -i $RANDBED -g $GENOME -s -l $FLANKSIZE -r 0 > $RANDBED".flank_sl1000r0.bed"

  

  
  
