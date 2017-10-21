#!/bin/bash

cat src/aligned_bam_to_fastq.cmds | parallel -j $1
