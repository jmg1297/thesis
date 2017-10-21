library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

run_ballgown <- function(metadata_file, data_dir, output_csv, p_max=0.05, q_max=0.05){
  metadata <- read.csv(metadata_file)
  bg_data <- ballgown(dataDir=data_dir, samplePattern="C57BL-6J", pData=metadata)
  # filtering step?
  bg_data_filt <- bg_data

  cell_de_results <- stattest(
    bg_data_filt,
    feature="transcript",
    covariate="cell",
    adjustvars=c("sex"),
    getFC=TRUE,
    meas="FPKM"
  )

  cell_de_results_transcripts <- data.frame(
    transcriptNames=ballgown::transcriptNames(bg_data_filt),
    transcriptIDs=ballgown::transcriptIDs(bg_data_filt),
    cell_de_results
  )

  cell_de_results_transcripts <- arrange(cell_de_results_transcripts, pval)

  cell_de_results_transcripts_sig <- subset(
    cell_de_results_transcripts,
    cell_de_results_transcripts$qval < q_max
  )

  write.csv(cell_de_results_transcripts_sig, output_csv)
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==4){
  run_ballgown(args[1], args[2], args[3], args[4])
} else if (length(args)==5){
  run_ballgown(args[1], args[2], args[3], p_max=args[4], q_max=args[5])
} else {
  run_ballgown(args[1], args[2], args[3])
}
