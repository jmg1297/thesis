#!/bin/bash

blastn \
-db $HOME/projects/cast_ref/genome/blast_db/CAST_EiJ.chromosomes.unplaced.gt2k.fa \
-query data/mm10_rrna_fasta/mm10_rRNA.fa \
-out data/blast_rRNA_to_cast/mm10_rRNA-CAST_genome.blast.out.plusonly \
-num_threads 4 \
-strand plus \
-outfmt 6
