#!/bin/bash

cat src/make_novel_beds.cmds | parallel -j $1
