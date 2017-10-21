#!/bin/bash

REF=$HOME/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf
for f in $(find data/kallisto_quant/ -type f -name "abundance.tsv");
do
SAMPLE=$(echo $f | cut -d/ -f3)
echo $SAMPLE
OUTDIR=data/ensembl_genes_kallisto_tpms/$SAMPLE
mkdir $OUTDIR
paste \
<(cat $REF | sed '/^#/d' | grep -P "\sgene\s" | awk '{print $1,$4,$5,$7,$10,$14,$18}' | sed -r 's/[;"]//g' | LC_ALL=C sort -k5 | sed '/^[GJM]/d') \
<(paste \
<(cat $REF | sed '/^#/d' | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$10,$14,$18,$22,$7}' | sed -r 's/[";]//g' | LC_ALL=C sort -k5 | sed '/^[GJM]/d') \
<(cat $f | sed '1d' | awk '{print $1,$5}' | LC_ALL=C sort -k1 | sed -r 's/_[0-9]+//' | sed -r 's/\s+/\t/g' | bedtools groupby -g 1 -c 2 -o sum) | \
awk '{print $1,$2,$3,$4,$6,$7,$10,$8}' | LC_ALL=C sort -k4 | sed -r 's/\s+/\t/g' | bedtools groupby -g 4 -c 7 -o sum | LC_ALL=C sort -k1) | \
awk '{print $1,$2,$3,$7,$5,$6,$9,$4}' | \
sed -r -e 's/\s+/:/4g' -e 's/:/ /3g' -e 's/\s+/\t/g' | \
sort -k1,1 -k2,2n \
> $OUTDIR/genes.bed
done

cd data/ensembl_genes_kallisto_tpms
cat $(find . -type f -name "genes.bed") | sort -k1,1 -k2,2n | bedtools groupby -g 1-4 -c 5 -o mean > genes_mean_tpms.bed
