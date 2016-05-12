#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
  cat('Usage: scriptName.R -p=output.pdf -file=file (header must be given, no row names)\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-column\tINT\t(Optional) Column to be used. Default: 1\n',file=stderr())
  cat('\t-breaks\tSTR\t(Optional) One of:\
                      \t\ta vector giving the breakpoints between histogram cells,
                      \t\ta function to compute the vector of breakpoints,
                      \t\ta single number giving the number of cells for the histogram,
                      \t\ta character string naming an algorithm to compute the number of cells,
                      \t\ta character string naming an algorithm to compute the number of cells. Default: Sturges\n', file=stderr())
  cat('\t-freq\tSTR\t(Optional) If true, the histogram graphic is a  representation of frequencies. Default: NULL\n', file=stderr())
  cat('\t-right\tSTR\t(Optional) If true, the histogram cells are right-closed (left open) intervals. Default: T\n', file=stderr())
  cat('\t-col\tSTR\t(Optional) A colour to be used to fill the bars. Default: NULL\n', file=stderr())
  cat('\t-p|-pdf\tFILE\tThe output figure in pdf [hist.pdf].\n',file=stderr())
  cat('\t-m|main\tSTR\t(Optional) The main title. Default: Main Title\n',file=stderr())
  cat('\t-xlab\tSTR\t(Optional) A title for x axis Default: X\n',file=stderr())
  cat('\t-ylab\tSTR\t(Optional) A title for y axis Default: Y\n',file=stderr())
  cat('\t-h\t\tShow help\n',file=stderr())
  q(save='no')
}

myColumn=1
myBreaks="Sturges"
myFreq=NULL
myRight='T'
myCol=NULL
myPdf='hist.pdf'
myMain='Main Title'
myXLab='X'
myYLab='Y'

if(length(args) >=1){
  for(i in 1:length(args)){
    arg=args[i]
    if(arg == '-h'){
      usage()
    }
    if(grepl('^-column?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -column')
      }else{
        myColumn=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-breaks?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -breaks')
      }else{
        myBreaks=arg.split[2]
      }
    }
    if(grepl('^-freq?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -freq')
      }else{
        myFreq=as.logical(arg.split[2])
      }
    }
    if(grepl('^-right?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -right')
      }else{
        myRight=as.logical(arg.split[2])
      }
    }
    if(grepl('^-col?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -col')
      }else{
        myCol=arg.split[2]
      }
    }
    if(grepl('^-p(df)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p')
      }else{
        myPdf=arg.split[2]
      }
    }
    if(grepl('^-m(ain)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -name2')
      }else{
        myMain=arg.split[2]
      }
    }
    if(grepl('^-xlab?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -xlab')
      }else{
        myXLab=arg.split[2]
      }
    }
    if(grepl('^-ylab?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -ylab')
      }else{
        myYLab=arg.split[2]
      }
    }
    if(grepl('^-file?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -file')
      }else{
        file=arg.split[2]
      }
    }
  }
}

data=read.table(file=file,header=TRUE, sep="\t")
data=data[,myColumn]
pdf(myPdf)
hist(data, breaks=myBreaks, freq=myFreq, right=myRight, col=myCol, main=myMain, xlab=myXLab, ylab=myYLab)

