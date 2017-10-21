#!/bin/bash

for f in $(find -type f -name "*_merge.gtf");
do
samples=$(echo $f | cut -d/ -f4 | sed -r 's/(.*)_merge.gtf/\1/')
echo "gffcompare -R -r /home/jg600/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.protein_coding.gtf -o data/gffcompare_ensembl/"$samples" "$f >> src/gffcompare_ensembl.cmds
done
