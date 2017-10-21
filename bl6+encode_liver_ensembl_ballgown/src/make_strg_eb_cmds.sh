#!/bin/bash

ENSEMBL=$HOME/projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf
SAMPLEDIR=$HOME/projects/rna-seq_transcriptome_reconstruction/data/rrna_filtered

for d in $(find $SAMPLEDIR -type d -name "*C57BL-6J*");
do
  SAMPLE=$(echo $d | cut -d/ -f8)
  mkdir data/stringtie_eb/$SAMPLE
  echo "stringtie -e -B -p 4 -G "$ENSEMBL" -o data/stringtie_eb/"$SAMPLE"/"$SAMPLE".gtf "$SAMPLEDIR"/"$SAMPLE"/rseqc_rrna.ex.sorted.bam" >> src/stringtie_eb.cmds
done
