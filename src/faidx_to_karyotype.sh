#!/bin/bash

cat CAST_EiJ.chromosomes.unplaced.gt2k.fa.fai | grep ^chr | awk '{print $1,$2}' | sed -r 's/chr//' | sort -k1n | awk '{if (NR==1){store=$_} else {print}} END {print store}' | awk '{print "chr","-","cast"$1,0,$2,"chr"$1}' > CAST_EiJ.chromosomes.unplaced.gt2k.karyotype.txt
