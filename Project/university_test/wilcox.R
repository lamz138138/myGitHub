#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
  cat('Usage: scriptName.R -file1=file1 -file2=file2\n',file=stderr())
  cat('Option:\n',file=stderr())
  cat('\t-column1\tINT\t(Optional) Column in file 1. Default: 1\n',file=stderr())
  cat('\t-column2\tINT\t(Optional) Column in file 2. Default: 1\n',file=stderr())
  cat('\t-file1\t\tFILE\tFile 1.\n',file=stderr())
  cat('\t-file2\t\tFILE\tFile 2.\n',file=stderr())
  cat('\t-myAlternative\tSTR\t(Optional) Alternative hypothesis: "two.sided", "greater" or "less". Default: two.sided\n',file=stderr())
  cat('\t-myMu\t\tINT\t(Optional) A number specifying an optional parameter used to form the null hypothesis. Default: 0\n',file=stderr())
  cat('\t-myPaired\tSTR\t(Optional) A logical indicating whether you want a paired test. Default: NULL\n',file=stderr())
  cat('\t-myExact\tSTR\t(Optional) A logical indicating whether an exact p-value should be computed. Default: FALSE\n',file=stderr())
  cat('\t-myCorrect\tSTR\t(Optional) A logical indicating whether an exact p-value should be computed. Default: FALSE\n',file=stderr())
  cat('\t-myConfInt\tSTR\t(Optional) A logical indicating whether a confidence interval should be computed. Default: FALSE\n',file=stderr())
  cat('\t-myConfLevel\tFloat\t(Optional) Confidence level of the interval. Default: 0.95\n',file=stderr())
  cat('\t-NA\t\tSTR\t(Optional) The NA string. Default: NA\n',file=stderr())
  cat('\t-h\t\t\tShow help\n',file=stderr())
  q(save='no')
}

column1=1
column2=1
myNa='NA'
myAlternative='two.sided'
myMu=0
myPaired='FALSE'
myExact=NULL
myCorrect='FALSE'
myConfInt='FALSE'
myConfLevel=0.95
if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    if(arg == '-h'){
      usage()
    }
    if(grepl('^-column1?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -column1')
      }else{
        column1=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-column2?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -column2')
      }else{
        column2=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-file1?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -file1')
      }else{
        file1=arg.split[2]
      }
    }
    if(grepl('^-file2?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -file2')
      }else{
        file2=arg.split[2]
      }
    }
    if(grepl('^-myAlternative?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myAlternative')
      }else{
        myAlternative=arg.split[2]
      }
    }
    if(grepl('^-myMu?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myMu')
      }else{
        myMu=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-myPaired?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myPaired')
      }else{
        myPaired=arg.split[2]
      }
    }
    if(grepl('^-myExact?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myExact')
      }else{
        myExact=arg.split[2]
      }
    }
    if(grepl('^-myCorrect?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myCorrect')
      }else{
        myCorrect=arg.split[2]
      }
    }
    if(grepl('^-myConfInt?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myConfInt')
      }else{
        myConfInt=arg.split[2]
      }
    }
    if(grepl('^-myConfLevel?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -myConfLevel')
      }else{
        myConfLevel=as.double(arg.split[2])
      }
    }
  }
}

data1=read.table(file=file1,header=FALSE, na.strings=myNa, sep="\t")
data1<-data1[, column1]
data2=read.table(file=file2,header=FALSE, na.strings=myNa, sep="\t")
data2<-data2[, column2]

if(length(data1)<2 || length(data2)<2){
  pValue="Data less than 2"
}else{
  pValue=wilcox.test(data1, data2, alternative=myAlternative, mu=myMu, paired=as.logical(myPaired), exact=myExact,
                     correct=as.logical(myCorrect),conf.int=as.logical(myConfInt), conf.level=myConfLevel)$p.value
}

cat(pValue)
cat('\n')
