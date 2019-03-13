#!/bin/env Rscript

library("optparse")
library("ggplot2")


option_list <- list(
  make_option("--columns", action="store", type="character", default="all",
              help="Column to be draw [default: %default].")
)

parser <- OptionParser(usage="%prog [options] input output", option_list=option_list,
                       description="This script was used to draw boxplot [header and rownames must be provided].")

args <- parse_args(parser, positional_arguments = 2)
opt <- args$options
input <- args$args[1]
output <- args$args[2]

data=read.delim(file=input, sep="\t", row.names=1)
Columns=colnames(data)
if(opt$columns!="all"){
  Columns=strsplit(opt$columns, ",", fixed = TRUE)[[1]]
}

drawData="NA"
myIndex=0
for(i in Columns){
  if(myIndex=="0"){
    drawData=cbind(data[,i], i)
    drawData=as.data.frame(drawData)
    colnames(drawData)=c("value", "type")
    myIndex=1
  }else{
    temp_data=cbind(data[,i], i)
    temp_data=as.data.frame(temp_data)
    colnames(temp_data)=c("value", "type")
    drawData=rbind(drawData, temp_data)
  }
}

drawData$value=as.numeric(drawData$value)
outputPlot = ggplot(drawData, aes(x=type, y=value)) + geom_boxplot()

pdf(output)
outputPlot
dev.off()


