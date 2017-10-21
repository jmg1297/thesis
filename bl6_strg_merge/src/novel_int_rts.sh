#!/bin/bash

cat src/novel_int_rts.cmds | parallel -j $1
