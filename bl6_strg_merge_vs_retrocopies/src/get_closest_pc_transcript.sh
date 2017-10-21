#!bin/bash

TRANSCRIPTS=/home/jg600/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.protein_coding_transcripts.bed
BASEDIR=$(pwd)
TARGETDIR=data/expr_retrocopies_vs_protein_coding_genes

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
  MRG=$(echo $d | cut -d/ -f3)
  mkdir $TARGETDIR/$MRG
  bedtools closest -d -a <(sort -k1,1 -k2,2n $d/expressed_retrocopies.bed) -b <(sort -k1,1 -k2,2n $TRANSCRIPTS) > $TARGETDIR/$MRG/expr_retrocopies-pc_transcripts.closest
  
  cd data/unexpr_retrocopies_with_rand_strand/$MRG
  for f in $(find . -type f -name "*.bed");
  do
    bedtools closest -d -a <(sort -k1,1 -k2,2n $f) -b <(sort -k1,1 -k2,2n $TRANSCRIPTS) > $f."pc_transcripts.closest"
  done
  cd $BASEDIR
done

RANDBED=data/expr_retrocopies_vs_protein_coding_genes/ALL/rand_rclike_regions.18000.bed

python $HOME/projects/src/get_random_regions.py \
$HOME/projects/mm10_ref/genome/mm10.chrom.sizes \
$HOME/projects/mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.bed \
50 \
18000 \
data/tmp.bed

sort -k1,1 -k2,2n data/tmp.bed > $RANDBED
rm data/tmp.bed

bedtools closest -d -a $RANDBED -b <(sort -k1,1 -k2,2n $TRANSCRIPTS) > $RANDBED".pc_transcripts.closest"
