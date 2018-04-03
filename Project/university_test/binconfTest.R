#!/bin/env Rscript

library("Hmisc")

args = commandArgs(trailingOnly=TRUE)

data=read.table(args[1],header=TRUE,row.names=1)

for(i in 1:nrow(data)){
  output_1=rownames(data)[i]
  for(j in 1:8){
    pvalue=binconf(data[i,j], data[i,j+8], 0.05, method="wilson")[1]
    output_1=cbind(output_1, pvalue)
  }
  print(output_1)
}


