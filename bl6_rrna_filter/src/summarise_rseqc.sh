#!/bin/bash

cat src/summarise_rseqc.cmds | parallel -j $1
