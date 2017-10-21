#!/bin/bash

INDEX=data/kallisto_idx/kallisto_strg_merge.idx
FQDIR=$HOME/thesis_projects/stringtie_transcriptome_rsem/data/rrna_filtered_fastq

mkdir data/kallisto_quant

for d in $(find $FQDIR -type d -name "*C57BL-6J*");
do
  SAMPLE=$(echo $d | cut -d/ -f8)
  echo $SAMPLE
  mkdir data/kallisto_quant/$SAMPLE
  echo "kallisto quant -i "$INDEX" -o data/kallisto_quant/"$SAMPLE" --fr-stranded -t 4 "$FQDIR"/"$SAMPLE"/sequence_1.fastq "$FQDIR"/"$SAMPLE"/sequence_2.fastq" >> src/kallisto_quant.cmds
done
