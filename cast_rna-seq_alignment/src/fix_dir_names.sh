#!/bin/bash

JSON=$HOME/projects/data/rna-seq.fastq/code2metadata.json

cd data/star_aligned

for D in $(find . -type d | sed '1d');
do
  echo $D
  CODE=$(echo $D | cut -d/ -f2)
  echo $CODE
  export CODE
  MD=$(cat $JSON | python -c "import os, sys, json; d=json.load(sys.stdin)[os.environ.get('CODE')]; out='_'.join([d['strain'], d['sex'], d['cell']]); print(out)")
  echo $MD
  NEWDIR=$D"_"$MD
  echo $NEWDIR
  mv $D $NEWDIR
done
