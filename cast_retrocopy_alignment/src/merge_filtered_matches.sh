#!/bin/bash

INPUT=$1
OUTPUT=$2

bedtools merge -i $INPUT -s -c 4 -o distinct > $OUTPUT
