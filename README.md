# Online Resources

## Introduction

This repository contains the Online Resources that supplement my PhD thesis,
"The Contribution of Retrotransposons to the Transcriptomes of Murine
Somatic Cells". These include scripts written to carry out analysis, and
resulting figures and statistics. Data files have been excluded, as they
represent the work of several collaborators as well as my own, and also
because of the large size of these files.

## Directory Structure

Each directory represents one analysis step, aside from a few specific
directories containing only data or scripts relevant to multiple analysis
steps (listed below). Each analysis directory contains a README file
describing the contents, and the following directories:
1. `data`: contains data files for this analysis step, usually organised
into subdirectories representing different stages in data processing
2. `src`: contains scripts used to carry out the analysis
3. `imgs`: contains images produced by scripts

The majority of Python scripts have been written using the `argparse` 
package, so running `python <script_name.py> -h` or
`python <script_name.py> --help` will display usage information. Scripts
written in other languages will have descriptive comments, unless
usage is self-explanatory.

## Special Directories
1. `data`: Contains raw RNA-seq FASTQ files, and data from other
sequencing experiments that has been analysed by collaborators.
2. `src`: Contains scripts relevant to multiple analysis steps, or
that are used to process reference files.
3. `mm10_ref`: Contains the mm10 reference genome and relevant 
annotations
4. `cast_ref`: Contains the CAST/Eij reference genome and relevant
annotations
