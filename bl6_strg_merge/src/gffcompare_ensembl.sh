#!/bin/bash

cat src/gffcompare_ensembl.cmds | parallel -j $1
