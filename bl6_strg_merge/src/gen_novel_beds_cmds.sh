#!/bin/bash

for d in $(find data/split_by_gffcmp_class -type d | sed '1d');
do
samples=$(echo $d | cut -d/ -f3)
mkdir data/novel_beds/$samples;
echo "cat "$d"/novel_non-coding.gtf | grep -P \"\stranscript\s\" | awk '{print \$1,\$4,\$5,\$12,1,\$18,\$7}' | sed -r 's/[\";]//g' | sed -r 's/\s+/\t/g' > data/novel_beds/"$samples"/novel_non-coding.bed" >> src/make_novel_beds.cmds
done
