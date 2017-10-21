#!/bin/bash

cat src/make_rsem_ref.cmds | parallel -j $1
