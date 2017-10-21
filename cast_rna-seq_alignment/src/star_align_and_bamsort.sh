#!/bin/bash


cat src/star_align_and_bamsort.cmds | parallel -j $1
