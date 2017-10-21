#!/bin/bash

cat src/kallisto_ensembl_quant.cmds | parallel -j $1
