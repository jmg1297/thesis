#!/bin/bash

cat src/stringtie.cmds | parallel -j $1
