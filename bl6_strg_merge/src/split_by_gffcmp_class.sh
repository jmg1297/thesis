#!/bin/bash

cat src/split_by_gffcmp_class.cmds | parallel -j $1
