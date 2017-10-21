#!/bin/bash

for d in $(find data/split_by_gffcmp_class -type d | sed '1d');
do
/usr/bin/python $HOME/projects/src/gtf_to_db.py $d/novel_non-coding.gtf
done
