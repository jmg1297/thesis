#!/bin/bash

GTF=$1
FASTA=$(echo $GTF | sed -r 's/gtf$/fa/')

echo "Extracting transcript sequences from "$GTF", writing to "$FASTA

GENOME=$2

gffread -w $FASTA -g $GENOME $GTF
