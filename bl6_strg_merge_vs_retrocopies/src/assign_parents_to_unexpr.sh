#!/bin/bash

BASEDIR=$(pwd)
DB=$HOME/projects/mm10_ref/retrogenes/retrogene_parent.db

for d in $(find data/all_exons_int_retrocopy_blocks -type d | sed '1d');
do
cd $d
/usr/bin/python $BASEDIR/src/append_parent_bed.py unexpressed_retrocopies.bed 4 $DB unexpressed_retrocopies_parents.mbed
cd $BASEDIR
done
