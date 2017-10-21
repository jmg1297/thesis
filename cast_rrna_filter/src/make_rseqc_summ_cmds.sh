#!/bin/bash

DIR=data/rrna_filter

for d in $(find $DIR -type d -name "*CAST*");
do
  echo "cd "$d"; echo 'hits,num_reads' > ex_hits_summary.csv; samtools view rseqc_rrna.ex.bam | awk '{print \$1,\$12}' | sort | uniq | awk '{print \$2}' | cut -d: -f3 | sort -n | uniq -c | awk '{print \$2,\$1}' | sed -r 's/\s+/,/g' >> ex_hits_summary.csv; samtools view rseqc_rrna.in.bam | awk '{print \$1}' | sort | uniq | wc -l > number_in.txt" >> src/summarise_rseqc.cmds
done
