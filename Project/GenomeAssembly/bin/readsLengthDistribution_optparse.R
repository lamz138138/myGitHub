#!/bin/env Rscript

library("optparse")

option_list <- list(
  make_option("--bins", action="store", type="integer", default=60,
              help="Number of bins [default: %default]. Overridden by binwidth"),
  make_option("--binWidth", action="store", type="double", 
              help="The width of the bins. Overridden by binwidth [default: null]"),
  make_option("--breaks", action="store", type="character", 
              help="Numeric vector giving the bin boundaries [default: null]"),
  make_option("--colName", action="store", type="character", default="len",
              help="The column name of length [default: %default]"),
  make_option("--vertical", action="store", type="double",
              help="The vertical line"),
  make_option("--xlimCoord", action="store", type="character",
              help="The range of x. Would add this to right top")
)

parser <- OptionParser(usage="%prog [options] input output", option_list=option_list,
                       description="This script was used to draw read length distribution.")

args <- parse_args(parser, positional_arguments = 2)
opt <- args$options
input <- args$args[1]
output <- args$args[2]

library("ggplot2")
library("grid")

data=read.delim(file=input, sep="\t", row.names=1)
outputPlot = ggplot(data)
if(!is.null(opt$breaks)){
	outputPlot = outputPlot + geom_histogram(aes(x=data[,opt$colName], fill=..count..),breaks=opt$breaks)
}else if(!is.null(opt$binWidth)){
	outputPlot = outputPlot + geom_histogram(aes(x=data[,opt$colName], fill=..count..),binwidth=opt$binWidth)
}else{
	outputPlot = outputPlot + geom_histogram(aes(x=data[,opt$colName], fill=..count..),bins=opt$bins) 
}

if(!is.null(opt$vertical)){
  outputPlot = outputPlot + geom_vline(xintercept=opt$vertical)
}

if(is.null(opt$xlimCoord)){
  pdf(output)
  outputPlot
  #dev.off()
}else{
  library("grid")
  pdf(output)
  print(outputPlot)
  outputPlot2 = outputPlot + coord_cartesian(xlim=as.numeric(strsplit(opt$xlimCoord,",")[[1]]))
  print(outputPlot2, vp=viewport(width = 0.4, height = 0.4, x=1, y=1, just = c("right", "top")))
  #dev.off()
}
