#!/bin/bash

DATADIR=data/intergenic
GFFCMPDIR=$HOME/thesis_projects/stringtie_merge/data/gffcompare_ensembl
ENSEMBLBED=$HOME/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.protein_coding_transcripts.bed
STRGRTDIR=$HOME/thesis_projects/strg_merge_all_exons_vs_rt_blocks/data/all_exons_int_rt_blocks
MBEDFN=annotated_u-ensembl_protein_coding.window+all_mstrg_matches+rt_overlaps.ext_mbed
BALLGOWN=$HOME/thesis_projects/ensembl_ballgown/data/ballgown/cell_de_results_transcripts.tsv
ENSEMBLTPM=$HOME/thesis_projects/kallisto_ensembl/data/kallisto_quant
STRGTPM=$HOME/thesis_projects/kallisto_strg_merge/data/kallisto_quant

for MRG in B T ALL;
do
  python src/get_class_code_gtf.py \
    $GFFCMPDIR/$MRG/$MRG.annotated.gtf \
    u \
    $DATADIR/$MRG.annotated_u.gtf

  cat $DATADIR/$MRG.annotated_u.gtf | \
    grep -P "\stranscript\s" | \
    awk '{print $1,$4,$5,$10,"0",$7}' | \
    sed -r -e 's/[\";]//g' -e 's/\s+/\t/g' | \
    sort -k1,1 -k2,2n \
    > $DATADIR/$MRG.annotated_u.bed

done

for WINDOW in 1000 5000 10000;
do
  for MRG in B T ALL;
  do
    bedtools window \
      -w $WINDOW \
      -a $DATADIR/$MRG.annotated_u.bed \
      -b $ENSEMBLBED \
      > $DATADIR/window_$WINDOW/$MRG.annotated_u-ensembl_protein_coding.window

    python src/append_matching_strg.py \
      $DATADIR/window_$WINDOW/$MRG.annotated_u-ensembl_protein_coding.window \
      $GFFCMPDIR/ALL/ALL.ALL_merge.gtf.tmap \
      $DATADIR/window_$WINDOW/$MRG.annotated_u-ensembl_protein_coding.window+all_mstrg_matches

    python src/append_rt_overlaps.py \
      $DATADIR/window_$WINDOW/$MRG.annotated_u-ensembl_protein_coding.window+all_mstrg_matches \
      $STRGRTDIR/$MRG/transcripts_retrotransposons.mbed \
      $STRGRTDIR/ALL/transcripts_retrotransposons.mbed \
      $DATADIR/window_$WINDOW/$MRG.$MBEDFN

    python src/split_mbed_by_rt_agreement.py \
      $DATADIR/window_$WINDOW/$MRG.$MBEDFN \
      $DATADIR/window_$WINDOW/$MRG.$MBEDFN

    mv $DATADIR/window_$WINDOW/$MRG.$MBEDFN $DATADIR/window_$WINDOW/$MRG.$MBEDFN".mbed"

    for STAT in "" no_rts. rt_agree. rt_disagree.;
    do
      cat $DATADIR/window_$WINDOW/$MRG.$MBEDFN.$STAT"mbed" | awk '{print $11}' | sort | uniq > $DATADIR/window_$WINDOW/$MRG.$STAT"ensembl.txt"
    done
  done
done

for WINDOW in 1000 5000 10000;
do
  for STAT in "" no_rts. rt_agree. rt_disagree.;
  do
    comm -12 $DATADIR/window_$WINDOW/B.$STAT"ensembl.txt" $DATADIR/window_$WINDOW/T.$STAT"ensembl.txt" > $DATADIR/window_$WINDOW/cell_shared.$STAT"ensembl.txt"
    comm -13 $DATADIR/window_$WINDOW/B.$STAT"ensembl.txt" $DATADIR/window_$WINDOW/T.$STAT"ensembl.txt" > $DATADIR/window_$WINDOW/T_spec.$STAT"ensembl.txt"
    comm -23 $DATADIR/window_$WINDOW/B.$STAT"ensembl.txt" $DATADIR/window_$WINDOW/T.$STAT"ensembl.txt" > $DATADIR/window_$WINDOW/B_spec.$STAT"ensembl.txt"

    python src/ensembl_fc_volcano.py \
      $DATADIR/window_$WINDOW/B_spec.$STAT"ensembl.txt" \
      $DATADIR/window_$WINDOW/T_spec.$STAT"ensembl.txt" \
      $DATADIR/window_$WINDOW/cell_shared.$STAT"ensembl.txt" \
      $BALLGOWN \
      "imgs/intergenic_w"$WINDOW"_"$STAT"ensembl_volcano.svg" \
      ""

    cat $DATADIR/window_$WINDOW/ALL.$MBEDFN.$STAT"mbed" | \
      awk '{print $4,$11}' | \
      sort | \
      uniq | \
      sed -r 's/\s+/,/g' \
      > $DATADIR/window_$WINDOW/ALL.$STAT"strg_ensembl_pairs.csv"

    python src/scatter_strg_mrg_vs_ensembl_tpm.py \
      $DATADIR/window_$WINDOW/ALL.$STAT"strg_ensembl_pairs.csv" \
      $ENSEMBLTPM \
      $STRGTPM \
      "imgs/intergenic_w"$WINDOW"_"$STAT".strg_ensembl_tpm_scatter.svg"

    python src/scatter_strg_mrg_vs_ensembl_tpm.py \
      -m 1 \
      $DATADIR/window_$WINDOW/ALL.$STAT"strg_ensembl_pairs.csv" \
      $ENSEMBLTPM \
      $STRGTPM \
      "imgs/intergenic_w"$WINDOW"_"$STAT"strg_ensembl_tpm_scatter.tpmgt1.svg"
  done
done
