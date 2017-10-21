#!/bin/bash

cat src/rt_content_hmap.cmds | parallel -j $1
