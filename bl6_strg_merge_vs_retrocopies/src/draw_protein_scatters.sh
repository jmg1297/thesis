#!/bin/bash

DATADIR=data/all_exons_int_retrocopy_blocks
GTFFNAME=/home/jg600/thesis_projects/mm10_ref/transcriptome/ensembl/Mus_musculus.GRCm38.84.gtf
ALLPROTTSV=/home/jg600/thesis_projects/data/proteomes/all/SPGRP_mouse_proteins.median_normed.tsv
FILTPROTTSV=/home/jg600/thesis_projects/data/proteomes/fdr_filtered/SPGRP_mouse_proteins_fdr-filtered.median_normed.tsv

# First scatter parents with cell-specific retrocopies
python src/scatter_cell_spec_parent_prot_levels.py \
$DATADIR/B_spec_parents_QUICKPATH.txt \
$DATADIR/T_spec_parents_QUICKPATH.txt \
$DATADIR/cell_shared_parents_QUICKPATH.txt \
$GTFFNAME $ALLPROTTSV imgs/cell_spec_parent_protein_scatter.svg

python src/scatter_cell_spec_parent_prot_levels.py \
$DATADIR/B_spec_parents_QUICKPATH.txt \
$DATADIR/T_spec_parents_QUICKPATH.txt \
$DATADIR/cell_shared_parents_QUICKPATH.txt \
$GTFFNAME $FILTPROTTSV imgs/cell_spec_parent_fdr-filt-protein_scatter.svg

# Now scatter those with upregulation in the presence of a retrocopy
python src/scatter_cell_spec_parent_prot_levels.py \
$DATADIR/B_spec_parents_with_B_upreg.txt \
$DATADIR/T_spec_parents_with_T_upreg.txt \
$DATADIR/cell_shared_parents_QUICKPATH.txt \
$GTFFNAME $ALLPROTTSV imgs/cell_spec+upreg_parent_protein_scatter.svg

python src/scatter_cell_spec_parent_prot_levels.py \
$DATADIR/B_spec_parents_with_B_upreg.txt \
$DATADIR/T_spec_parents_with_T_upreg.txt \
$DATADIR/cell_shared_parents_QUICKPATH.txt \
$GTFFNAME $FILTPROTTSV imgs/cell_spec+upreg_parent_fdr-filt-protein_scatter.svg
