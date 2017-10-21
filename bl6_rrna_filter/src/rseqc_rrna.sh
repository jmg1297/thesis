#!/bin/bash

cat src/rseqc_rrna.cmds | parallel -j $1
