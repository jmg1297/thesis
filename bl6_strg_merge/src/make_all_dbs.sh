#!/bin/bash

for f in $(find data/stringtie_merge -type f -name "*_merge.gtf");
do
/usr/bin/python $HOME/projects/src/gtf_to_db.py $f
done
