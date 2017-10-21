#!/bin/bash
# One liner to create a BED file of Ensembl genes with their ID, name, and biotype, along with the TPM for that gene from Kallisto.

# Steps:
# Paste together:
# 	Relevant fields for genes from the Ensembl annotation GTF, sorted by gene ID
#	Paste together:
#		Relevant fields for transcripts from the Ensembl annotation GTF, sorted by transcript ID
#		Kallisto TPMs, deduplicated, collapsed, and sorted by transcript ID
#	then extract relevant fields, sort by gene ID, collapse into genes and sum TPMs, sort by gene ID
# then extract relevant fields, format to BED specification and sort by coordinate
# 

paste \
<(cat Mus_musculus.GRCm38.84.gtf | sed '/^#/d' | grep -P "\sgene\s" | awk '{print $1,$4,$5,$7,$10,$14,$18}' | sed -r 's/[;"]//g' | LC_ALL=C sort -k5 | sed '/^[GJM]/d') \
<(paste \
<(cat Mus_musculus.GRCm38.84.gtf | sed '/^#/d' | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$10,$14,$18,$22,$7}' | sed -r 's/[";]//g' | LC_ALL=C sort -k5 | sed '/^[GJM]/d') \
<(cat abundance.tsv | sed '1d' | awk '{print $1,$5}' | LC_ALL=C sort -k1 | sed -r 's/_[0-9]+//' | sed -r 's/\s+/\t/g' | bedtools groupby -g 1 -c 2 -o sum) | \
awk '{print $1,$2,$3,$4,$6,$7,$10,$8}' | LC_ALL=C sort -k4 | sed -r 's/\s+/\t/g' | bedtools groupby -g 4 -c 7 -o sum | LC_ALL=C sort -k1) | \
awk '{print $1,$2,$3,$7,$5,$6,$9,$4}' | \
sed -r -e 's/\s+/:/4g' -e 's/:/ /3g' -e 's/\s+/\t/g' | \
sort -k1,1 -k2,2n | \
less -S 

