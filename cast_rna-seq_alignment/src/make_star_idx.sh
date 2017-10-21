#!/bin/bash

STAR \
--runThreadN $1 \
--runMode genomeGenerate \
--genomeDir data/STAR_index \
--genomeFastaFiles $HOME/projects/cast_ref/genome/blast_db/CAST_EiJ.chromosomes.unplaced.gt2k.fa 
