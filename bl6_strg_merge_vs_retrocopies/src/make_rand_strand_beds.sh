#!/bin/bash

python src/assign_random_strands.py data/all_exons_int_retrocopy_blocks/ALL/unexpressed_retrocopies_parents.mbed data/all_exons_int_retrocopy_blocks/ALL/strand_combo_counts.json test.bed

TARGETDIR=data/unexpr_retrocopies_with_rand_strand

for d in $(find data/all_exons_int_retrocopy_blocks/ -type d | sed '1d');
do
MRG=$(echo $d | cut -d/ -f3)
mkdir $TARGETDIR/$MRG
for i in $(seq 1 $1);
do
python src/assign_random_strands.py $d/unexpressed_retrocopies_parents.mbed $d/strand_combo_counts.json $TARGETDIR/$MRG/unexpressed_retrocopies_rand_strand.set$i.bed
done
done

