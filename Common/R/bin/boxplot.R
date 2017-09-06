#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
  cat('Usage: scriptName.R -p=output.pdf -file=file (header must be given, no row names)\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-column\t\tINT\t(Optional) Column to be used. Default: 1\n',file=stderr())
  cat('\t-col\t\tSTR\t(Optional) A colour to be used to fill the bars. Default: NULL\n', file=stderr())
  cat('\t-p|-pdf\t\tFILE\tThe output figure in pdf [hist.pdf].\n',file=stderr())
  cat('\t-m|main\t\tSTR\t(Optional) The main title. Default: Main Title\n',file=stderr())
  cat('\t-xlab\t\tSTR\t(Optional) A title for x axis Default: X\n',file=stderr())
  cat('\t-ylab\t\tSTR\t(Optional) A title for y axis Default: Y\n',file=stderr()) 
  cat('\t-outline\tSTR\t(Optional) Whether draw outliers. Default: FALSE\n', file=stderr())
  cat('\t-names\t\tSTR\t(Optional) Group labels which will be printed under each boxplot\n', file=stderr())
  cat('\t-log\t\tSTR\t(Optional) if x or y or both coordinates should be plotted in log scale, can be "x", "y" or ?\n', file=stderr())
  cat('\t-horizontal\tSTR\t(Optional) If the boxplots should be horizontal. Default: FALSE\n', file=stderr())
  cat('\t-h\t\t\tShow help\n',file=stderr())
  q(save='no')
}

myColumn=1
myCol=NULL
myPdf='boxplot.pdf'
myMain='Main Title'
myXLab='X'
myYLab='Y'
myOutline=as.logical('F')
myNames='name'
myLog='x'
myHorizontal=as.logical('F')

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
    if(grepl('^-outline?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -outline')
      }else{
        myOutline=as.logical(arg.split[2])
      }
    }
    if(grepl('^-names?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -names')
      }else{
        myNames=arg.split[2]
      }
    }
    if(grepl('^-log?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -log')
      }else{
        myLog=arg.split[2]
      }
    }
    if(grepl('^-horizontal?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -horizontal')
      }else{
        myHorizontal=as.logical(arg.split[2])
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
boxplot(data, outline=myOutline, horizontal=myHorizontal, main=myMain, xlab=myXLab, ylab=myYLab)