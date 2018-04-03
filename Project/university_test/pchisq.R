#!/bin/env Rscript

library("metap")

args = commandArgs(trailingOnly=TRUE)

data=read.table(args[1],header=FALSE,sep="\t")

for(i in 1:nrow(data)){
  tranName=as.character(data[i,1])
  geneName=as.character(data[i,2])
  myValue=pchisq((sum(log(data[i,3:34]))*-2),df=64,lower.tail=F)
  cat(tranName, geneName, as.character(myValue), sep="\t", fill=TRUE)
}
