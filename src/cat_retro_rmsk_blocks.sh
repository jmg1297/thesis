#!/bin/bash

cat mm10_ref/repeatmasker/formatted_beds/rmskJoinedBaseline.blocks.formatted.bed mm10_ref/retrogenes/formatted_beds/ucscRetroInfo6.blocks.bed | \
sort -k1,1 -k2,2n \
> mm10_ref/rmskJoinedBaseline_blocks+ucscsRetroInfo6_blocks.formatted.bed
