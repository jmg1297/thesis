#!/bin/bash

makeblastdb -in cast_ref/genome/CAST_EiJ.chromosomes.unplaced.gt2k.fa -parse_seqids -dbtype nucl

blastn \
-db cast_ref/genome/blast_db/CAST_EiJ.chromosomes.unplaced.gt2k.fa \
-query mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.unstranded.blocks.fa \
-out cast_ref/genome/mm10_ucscRetroInfo6-CAST_genome.blast.out \
-num_threads 8 \
-outfmt 6
