#!/bin/bash

cat src/rsem_expression.cmds | parallel -j $1
