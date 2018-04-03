#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

data=read.table(args[1],header=FALSE,sep="\t")

adjustData=p.adjust(data[,3],method="BH", n=nrow(data))

for(i in 1:nrow(data)){
  cat(as.character(data[i,1]), as.character(data[i,2]), data[i,3], adjustData[i], sep="\t", fill=TRUE)
}
