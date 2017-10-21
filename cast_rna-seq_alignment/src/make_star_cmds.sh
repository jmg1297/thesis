#!/bin/bash

SAMPLES=$HOME/projects/data/rna-seq.fastq/CAST_samples.txt
DATADIR=$HOME/projects/data/rna-seq.fastq
GENOMEDIR=data/STAR_index

for s in $(cat $SAMPLES);
do
LANE=$(echo $s | cut -d_ -f2 | cut -d# -f1)
FQ1=$DATADIR"/"$s"/s_"$LANE"_1_sequence.fq"
FQ2=$DATADIR"/"$s"/s_"$LANE"_2_sequence.fq"
OUTPREFIX=data/star_aligned/$s/
mkdir $OUTPREFIX

echo "STAR --genomeDir "$GENOMEDIR" --runThreadN 4 --outReadsUnmapped Fastx --outSAMattributes All --outFilterMultimapNmax 50 --readFilesIn "$FQ1" "$FQ2" --outFileNamePrefix "$OUTPREFIX"; cd "$OUTPREFIX"; samtools view -@ 4 -Sb Aligned.out.sam > Aligned.out.bam; samtools sort -@ 4 -T tmp.sorted -O 'bam' Aligned.out.bam > Aligned.out.sorted.bam;" >> src/star_align_and_bamsort.cmds

done
