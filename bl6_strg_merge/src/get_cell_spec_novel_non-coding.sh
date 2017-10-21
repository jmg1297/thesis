#!/bin/bash

mkdir data/cell_spec_novel_non-coding/refB_queryT
gffcompare -R -r data/split_by_gffcmp_class/B/novel_non-coding.gtf -o data/cell_spec_novel_non-coding/refB_queryT/rBqT data/split_by_gffcmp_class/T/novel_non-coding.gtf
mv data/split_by_gffcmp_class/T/*gtf.* data/cell_spec_novel_non-coding/refB_queryT/

mkdir data/cell_spec_novel_non-coding/refT_queryB
gffcompare -R -r data/split_by_gffcmp_class/T/novel_non-coding.gtf -o data/cell_spec_novel_non-coding/refT_queryB/rTqB data/split_by_gffcmp_class/B/novel_non-coding.gtf
mv data/split_by_gffcmp_class/B/*gtf.* data/cell_spec_novel_non-coding/refT_queryB/

#----------------------------------------------------------------------------------------------------------------------------------

mkdir data/cell_spec_novel_non-coding/refB_queryT/split_by_class
python src/split_by_gffcompare.py \
data/cell_spec_novel_non-coding/refB_queryT/rBqT.novel_non-coding.gtf.tmap \
data/split_by_gffcmp_class/T/novel_non-coding.gtf \
data/cell_spec_novel_non-coding/refB_queryT/split_by_class

cat data/cell_spec_novel_non-coding/refB_queryT/split_by_class/novel_non-coding.gtf | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$12,1,$18,$7}' | sed -r 's/[";]//g' | sed -r 's/\s+/\t/g' \
> data/cell_spec_novel_non-coding/refB_queryT/split_by_class/novel_non-coding.bed

mkdir data/cell_spec_novel_non-coding/refT_queryB/split_by_class
python src/split_by_gffcompare.py \
data/cell_spec_novel_non-coding/refT_queryB/rTqB.novel_non-coding.gtf.tmap \
data/split_by_gffcmp_class/B/novel_non-coding.gtf \
data/cell_spec_novel_non-coding/refT_queryB/split_by_class

cat data/cell_spec_novel_non-coding/refT_queryB/split_by_class/novel_non-coding.gtf | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$12,1,$18,$7}' | sed -r 's/[";]//g' | sed -r 's/\s+/\t/g' \
> data/cell_spec_novel_non-coding/refT_queryB/split_by_class/novel_non-coding.bed

#------------------------------------------------------------------------------------------------------------------------------------

mkdir data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons
bedtools intersect -wo \
-a data/cell_spec_novel_non-coding/refB_queryT/split_by_class/novel_non-coding.bed \
-b /home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.formatted.bed \
> data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect

python $HOME/projects/stringtie_novel_vs_retrotransposons/src/transcript_intersect_hmap.py \
--intersectFilename data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect \
--outputImg data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/rt_overlap_hmap.svg \
--sortBy cluster


mkdir data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons
bedtools intersect -wo \
-a data/cell_spec_novel_non-coding/refT_queryB/split_by_class/novel_non-coding.bed \
-b /home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.formatted.bed \
> data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect

python $HOME/projects/stringtie_novel_vs_retrotransposons/src/transcript_intersect_hmap.py \
--intersectFilename data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/novel_non-coding.rmskJoinedBaseline_RTs_formatted.intersect \
--outputImg data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/rt_overlap_hmap.svg \
--sortBy cluster

#------------------------------------------------------------------------------------------------------------------------------------

cat data/cell_spec_novel_non-coding/refB_queryT/split_by_class/annotated.gtf | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$12,1,$18,$7}' | sed -r 's/[";]//g' | sed -r 's/\s+/\t/g' \
> data/cell_spec_novel_non-coding/refB_queryT/split_by_class/annotated.bed

bedtools intersect -wo \
-a data/cell_spec_novel_non-coding/refB_queryT/split_by_class/annotated.bed \
-b /home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.formatted.bed \
> data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/annotated.rmskJoinedBaseline_RTs_formatted.intersect

python $HOME/projects/stringtie_novel_vs_retrotransposons/src/transcript_intersect_hmap.py \
--intersectFilename data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/annotated.rmskJoinedBaseline_RTs_formatted.intersect \
--outputImg data/cell_spec_novel_non-coding/refB_queryT/novel_int_retrotransposons/cell_shared_rt_overlap_hmap.svg \
--sortBy cluster

cat data/cell_spec_novel_non-coding/refT_queryB/split_by_class/annotated.gtf | grep -P "\stranscript\s" | awk '{print $1,$4,$5,$12,1,$18,$7}' | sed -r 's/[";]//g' | sed -r 's/\s+/\t/g' \
> data/cell_spec_novel_non-coding/refT_queryB/split_by_class/annotated.bed

bedtools intersect -wo \
-a data/cell_spec_novel_non-coding/refT_queryB/split_by_class/annotated.bed \
-b /home/jg600/projects/mm10_ref/repeatmasker/retrotransposons_only/formatted_beds/rmskJoinedBaseline.RTs.formatted.bed \
> data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/annotated.rmskJoinedBaseline_RTs_formatted.intersect


python $HOME/projects/stringtie_novel_vs_retrotransposons/src/transcript_intersect_hmap.py \
--intersectFilename data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/annotated.rmskJoinedBaseline_RTs_formatted.intersect \
--outputImg data/cell_spec_novel_non-coding/refT_queryB/novel_int_retrotransposons/cell_shared_rt_overlap_hmap.svg \
--sortBy cluster
