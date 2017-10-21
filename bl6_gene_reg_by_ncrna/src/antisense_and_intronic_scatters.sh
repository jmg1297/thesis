#!/bin/bash

ENSEMBLTPM=$HOME/thesis_projects/kallisto_ensembl/data/kallisto_quant
STRGTPM=$HOME/thesis_projects/kallisto_strg_merge/data/kallisto_quant

python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_intronic/intronic_pairs.csv $ENSEMBLTPM $STRGTPM imgs/intronic_tpm_scatter.svg
python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_intronic/intronic_pairs.csv $ENSEMBLTPM $STRGTPM imgs/intronic_tpm_scatter.tpmgt1.svg -m 1

python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_rt_intronic/with_rt.intronic_pairs.csv $ENSEMBLTPM $STRGTPM imgs/with_rt_intronic_tpm_scatter.svg
python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_rt_intronic/with_rt.intronic_pairs.csv $ENSEMBLTPM $STRGTPM imgs/with_rt_intronic_tpm_scatter.tpmgt1.svg -m 1

python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_antisense/antisense_pairs.csv $ENSEMBLTPM $STRGTPM imgs/antisense_tpm_scatter.svg
python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_antisense/antisense_pairs.csv $ENSEMBLTPM $STRGTPM imgs/antisense_tpm_scatter.tpmgt1.svg -m 1

python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_rt_antisense/with_rt.antisense_pairs.csv $ENSEMBLTPM $STRGTPM imgs/with_rt_antisense_scatter.svg
python src/scatter_strg_mrg_vs_ensembl_tpm.py data/with_rt_antisense/with_rt.antisense_pairs.csv $ENSEMBLTPM $STRGTPM imgs/with_rt_antisense_scatter.tpmgt1.svg -m 1
