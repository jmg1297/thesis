#!/bin/bash

RRNA_BED=$HOME/projects/mm10_ref/transcriptome/rRNA/mm10_rRNA.chr.bed
GENOME=$HOME/projects/mm10_ref/genome/mm10_genome.fa

bedtools getfasta -name -fo >(fold -w 60 > data/mm10_rrna_fasta/mm10_rRNA.fa) -fi $GENOME -bed $RRNA_BED

