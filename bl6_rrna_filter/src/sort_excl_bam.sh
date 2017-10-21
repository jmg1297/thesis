#!/bin/bash

cat src/sort_excl_bam.cmds | parallel -j $1
