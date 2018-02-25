#!/usr/bin/env/Rscript

print_matrix <- function(m) {
  cat("MATRIX:")
  cat(deparse(substitute(m)), sep="\n")
  for (s in apply(m, 1, paste, collapse=",")){
    cat(s, sep="\n")
  }
}

suppressPackageStartupMessages(library("sigclust2"))

args <- commandArgs(trailingOnly=TRUE)

data_csv <- args[1]

df <- read.csv(data_csv, header=FALSE)
dat <- as.matrix(df, byrow=TRUE)

shc_result <- shc(dat, metric="euclidean", linkage="ward.D2")

merge <- shc_result$hc_dat$merge
height <- as.matrix(shc_result$hc_dat$height)
order <- as.matrix(shc_result$hc_dat$order)
p_vals <- as.matrix(shc_result$p_norm)

print_matrix(merge)
print_matrix(height)
print_matrix(order)
print_matrix(p_vals)

# Need to convert this into a format Python can process and return it
# Could simply stream it line by line and use subprocess.check_output