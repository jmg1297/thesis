#!/bin/bash

for d in $(find data/gffcompare_ensembl -type d | sed '1d');
do
samples=$(echo $d | cut -d/ -f3)
gtf_fname=$samples"_merge.gtf"
mkdir data/split_by_gffcmp_class/$samples
echo "python src/split_by_gffcompare.py "$d"/"$samples"."$gtf_fname".tmap data/stringtie_merge/"$gtf_fname" data/split_by_gffcmp_class/"$samples >> src/split_by_gffcmp_class.cmds
done
