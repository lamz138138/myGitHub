#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
#	make_option("--outputDir", action="store", type="character",
#              help="Directory of output.")
)

parser <- OptionParser(usage="cat file | %prog", option_list=option_list, description="This script was used to output min, max, median and mean value, and read from stdin.")
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
#input <- args$args[1]
#output <- args$args[2]

d<-scan("stdin", quiet=TRUE)
cat(paste("min:",min(d),sep=" "), paste("max:",max(d),sep=" "), paste("median:",median(d),sep=" "), paste("mean:",mean(d),sep=" "), sep="\n")
