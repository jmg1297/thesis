#!/bin/bash

cd data/strg_merge_int_rt_indels

comm -12 <(sort B/expressed_rtis.bed) <(sort T/expressed_rtis.bed) | sort -k1,1 -k2,2n > expr_rtis.cell_shared.bed

comm -13 <(sort B/expressed_rtis.bed) <(sort T/expressed_rtis.bed) | sort -k1,1 -k2,2n > expr_rtis.T_spec.bed

comm -23 <(sort B/expressed_rtis.bed) <(sort T/expressed_rtis.bed) | sort -k1,1 -k2,2n > expr_rtis.B_spec.bed


