#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
	cat('Usage: scriptName.R -p=output.pdf -file=file (header and row names must be given)\n',file=stderr())
	cat('Option:\n',file=stderr())
	cat('\t-file\t\t\t\tFILE\tFile.\n',file=stderr())
	cat('\t-p|-pdf\t\t\t\tFILE\tThe output name of file. Default: heatmap2.pdf].\n',file=stderr())
    cat('\t-m|main\t\t\t\tSTR\t(Optional) The main title. Default: Main Title\n',file=stderr())
    cat('\t-NA\t\t\t\tSTR\t(Optional) The NA string. Default: NA\n',file=stderr())
    cat('\t-naRM\t\t\t\tSTR\t(Optional) Whether NAs should be removed. Default: FALSE\n',file=stderr())
    cat('\t-scale\t\t\t\tSTR\t(Optional) Whether value should be centered, can be "row" or "column" or "none". Default: none\n',file=stderr())
	cat('\t-Rowv\t\t\t\tSTR\t(Optional) If rows should be clustered, can be "NULL", "FALSE" or "TRUE".  Default: TRUE\n',file=stderr())
	cat('\t-Colv\t\t\t\tSTR\t(Optional) If columns should be clustered, can be "NULL", "FALSE", "TRUE" or "Rowv".  Default: TRUE\n',file=stderr())
	cat('\t-distfun\t\t\tSTR\t(Optional) Function used to compute the distance (dissimilarity) between both rows and columns. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".  Default: euclidean\n',file=stderr())
	cat('\t-hclustfun\t\t\tSTR\t(Optional) Function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. This must be one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
                                \t\t\t"median" (= WPGMC) or "centroid" (= UPGMC). Default: complete\n',file=stderr())
	cat('\t-dendrogram\t\t\tSTR\t(Optional) Character string indicating whether to draw "none", "row", "column" or "both" dendrograms. If Rowv (or Colv) is FALSE or NULL and dendrogram is "both", then a warning is issued and Rowv (or Colv) arguments are
                                \t\t\thonoured. Defaults: "both".\n',file=stderr())
    cat('\t-distfun2\t\t\tSTR\t(Optional) Whether calculated distance for column and row respectively, can be "TRUE" or "FALSE". Defalut: FALSE\n',file=stderr())
    cat('\t-distfun2Row\t\t\tSTR\t(Optional) Distance measure used in clustering rows when distfun2 is "TRUE", can be "pearson", "kendall" or "spearman". Default: pearson\n',file=stderr())
    cat('\t-distfun2Col\t\t\tSTR\t(Optional) Distance measure used in clustering columns when distfun2 is "TRUE", can be "pearson", "kendall" or "spearman". Default: pearson\n',file=stderr())
    cat('\t-hclustfun2Row\t\t\tSTR\t(Optional) The agglomeration method to be used in row when when distfun2 is "TRUE". This must be one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
                                \t\t\t"median" (= WPGMC) or "centroid" (= UPGMC). Default: complete\n',file=stderr())
    cat('\t-hclustfun2Col\t\t\tSTR\t(Optional) The agglomeration method to be used in column when when distfun2 is "TRUE". This must be one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
                                \t\t\t"median" (= WPGMC) or "centroid" (= UPGMC). Default: complete\n',file=stderr())                            
	cat('\t-width\t\t\t\tINT\t(Optional) The width in inches of output. Default: 1024\n',file=stderr())
	cat('\t-height\t\t\t\tINT\t(Optional) The height in inches of output. Default: 768\n',file=stderr())
	cat('\t-h\t\t\t\t\tShow help\n',file=stderr())
	q(save='no')
}

myPdf='heatmap2.pdf'
myNa='NA'
myMain='Main Title'
myNaRm='FALSE'
myScale='none'
myRowv='TRUE'
myColv='TRUE'
myDistfun='euclidean'
myHclustfun='complete'
myDendrogram='both'
myDistfun2='FALSE'
myDistfun2Row='pearson'
myDistfun2Col='pearson'
myHclustfun2Row='complete'
myHclustfun2Col='complete'
myWidth=7
myHeight=7
if(length(args) >= 1){
  for(i in 1:length(args)){
    arg=args[i]
    if(arg == '-h'){
      usage()
    }
    if(grepl('^-p(df)?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -p')
      }else{
        myPdf=arg.split[2]
      }
    }
    if(grepl('^-NA?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -NA')
      }else{
        myNa=arg.split[2]
      }
    }
    if(grepl('^-naRM?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -naRM')
      }else{
        myNaRm=arg.split[2]
      }
    }
    if(grepl('^-scale?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -scale')
      }else{
        myScale=arg.split[2]
      }
    }
    if(grepl('^-Rowv?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -Rowv')
      }else{
        myRowv=arg.split[2]
      }
    }
    if(grepl('^-Colv?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -Colv')
      }else{
        myColv=arg.split[2]
      }
    }
    if(grepl('^-distfun?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -distfun')
      }else{
        myDistfun=arg.split[2]
      }
    }
    if(grepl('^-hclustfun?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -hclustfun')
      }else{
        myHclustfun=arg.split[2]
      }
    }
    if(grepl('^-dendrogram?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -dendrogram')
      }else{
        myDendrogram=arg.split[2]
      }
    }
    if(grepl('^-width?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -width')
      }else{
        myWidth=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-height?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -height')
      }else{
        myHeight=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-distfun2?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -distfun2')
      }else{
        myDistfun2=arg.split[2]
      }
    }
    if(grepl('^-distfun2Row?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -distfun2Row')
      }else{
        myDistfun2Row=arg.split[2]
      }
    }
    if(grepl('^-distfun2Col?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -distfun2Col')
      }else{
        myDistfun2Col=arg.split[2]
      }
    }
    if(grepl('^-hclustfun2Row?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -hclustfun2Row')
      }else{
        myHclustfun2Row=arg.split[2]
      }
    }
    if(grepl('^-hclustfun2Col?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -hclustfun2Col')
      }else{
        myHclustfun2Col=arg.split[2]
      }
    }
    if(grepl('^-file?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -file')
      }else{
        myFile=arg.split[2]
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
  }
}

library(gplots)
data=read.table(file=myFile,header=TRUE, sep="\t",na.strings=myNa, row.names=1)
data=as.matrix(data)
if(myDistfun2=='FALSE'){
    distfunFunction<-function(x){
        dist(x, method=myDistfun)
    }
    hclustfunFuction<-function(x){
        hclust(x, method=myHclustfun)
    }
    pdf(myPdf, width=myWidth, height=myHeight)
    heatmap.2(data, na.rm=as.logical(myNaRm), scale=myScale, Rowv=as.logical(myRowv), Colv=as.logical(myColv), distfun=distfunFunction, hclustfun=hclustfunFuction, dendrogram=myDendrogram)
}else{
    rowCluster <- hclust(as.dist(1-cor(t(data), method=myDistfun2Row)), method=myHclustfun2Row)
    columnCluster <- hclust(as.dist(1-cor(data, method=myDistfun2Col)), method=myHclustfun2Col)
    pdf(myPdf, width=myWidth, height=myHeight)
    heatmap.2(data, na.rm=as.logical(myNaRm), scale=myScale, Rowv=as.dendrogram(rowCluster), Colv=as.dendrogram(columnCluster), density.info="none", trace="none")
}

