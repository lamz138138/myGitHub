#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
  cat('This script was used to draw read length distribution.\n', file=stderr())
  cat('Usage: scriptName.R -p=output.pdf -file=file\n', file=stderr())
  cat('Options:\n', file=stderr())
  cat('\t-file\t\tFILE\tFile (First row should be row name, and should contain header).\n', file=stderr())
  cat('\t-p|-pdf\t\tFILE\tThe output name of file. Default: readsLengthDis.pdf.\n', file=stderr())
  cat('\t-b|-bins\tINT\tNumber of bins. Overridden by binwidth. Default: 60.\n', file=stderr())
  cat('\t-binWidth\tINT\tThe width of the bins. Overridden by binwidth. Default: null.\n', file=stderr())
  cat('\t-breaks\tSTR\tNumeric vector giving the bin boundaries. Default: null.\n', file=stderr())
  cat('\t-colName\tSTRING\tThe column name of length. Default: len.\n', file=stderr())
  cat('\t-width\t\tINT\t(Optional) The width in inches of output. Default: 1024.\n',file=stderr())
  cat('\t-height\t\tINT\t(Optional) The height in inches of output. Default: 768.\n',file=stderr())
  cat('\t-m|main\t\tSTR\t(Optional) The main title. Default: Main Title.\n',file=stderr())
  cat('\t-xlab\t\tSTR\t(Optional) A title for x axis Default: X.\n',file=stderr())
  cat('\t-ylab\t\tSTR\t(Optional) A title for y axis Default: Y.\n',file=stderr())
  cat('\t-h\t\t\tShow help\n',file=stderr())
  q(save='no')
}

myPdf='readsLengthDis.pdf'
myWidth=1024
myHeight=768
myBins=60
myColName='len'
myMain='Main Title'
myXLab='X'
myYLab='Y'
myBinWidth='null'
myBreaks='null'

if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    # help
    if(arg == '-h'){
      usage()
    }
    
    # name of output
    if(grepl('^-p(df)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p')
      }else{
        myPdf=arg.split[2]
      }
    }
    
    # number of bins
    if(grepl('^-b(ins)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -b')
      }else{
        myBins=as.numeric(arg.split[2])
      }
    }
    
    # width of bin
    if(grepl('^-binWidth?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -binWidth')
      }else{
        myBinWidth=as.numeric(arg.split[2])
      }
    }

    # breaks
    if(grepl('^-breaks?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -breaks')
      }else{
        myBreaks=as.numeric(arg.split[2])
      }
    }

    # the column name of length
    if(grepl('^-colName?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -colName')
      }else{
        myColName=arg.split[2]
      }
    }
    
    # width of picture
    if(grepl('^-width?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -width')
      }else{
        myWidth=as.numeric(arg.split[2])
      }
    }
    
    # height of picture
    if(grepl('^-height?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -height')
      }else{
        myHeight=as.numeric(arg.split[2])
      }
    }
    
    # input file
    if(grepl('^-file?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -file')
      }else{
        myFile=arg.split[2]
      }
    }
    
    # main title
    if(grepl('^-m(ain)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -main')
      }else{
        myMain=arg.split[2]
      }
    }
    
    # X label
    if(grepl('^-xlab?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -xlab')
      }else{
        myXLab=arg.split[2]
      }
    }
    
    # Y label
    if(grepl('^-ylab?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -ylab')
      }else{
        myYLab=arg.split[2]
      }
    }
    
  }  
}

library("ggplot2")
data=read.delim(file=myFile, sep="\t", row.names=1)
#data=data[-length(data[,1]),]
outputPlot = ggplot(data)
if(myBreaks!='null'){
	outputPlot = outputPlot + geom_histogram(aes(x=data[,myColName]),breaks=myBreaks)
}else if(myBinWidth!='null'){
	outputPlot = outputPlot + geom_histogram(aes(x=data[,myColName]),binwidth=myBinWidth)
}else{
	outputPlot = outputPlot + geom_histogram(aes(x=data[,myColName]),bins=myBins) 
}
outputPlot = outputPlot + labs(title=myMain, x=myXLab, y=myYLab)
#ggsave(myPdf, width=myWidth, height=myHeight, plot=outputPlot, limitsize=FALSE)
ggsave(myPdf, plot=outputPlot)
