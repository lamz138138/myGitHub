#!/bin/env Rscript

args <- commandArgs(TRUE)
usage = function(){
	cat('Usage: scriptName.R -p=output.pdf -file=file (header and row names must be given)\n',file=stderr())
	cat('Option:\n',file=stderr())
	cat('\t-file\t\t\t\tFILE\tFile.\n',file=stderr())
	cat('\t-p|-pdf\t\t\t\tFILE\tThe output figure in pdf[pheatmap.pdf].\n',file=stderr())
	cat('\t-m|main\t\t\t\tSTR\t(Optional) The main title. Default: Main Title\n',file=stderr())
	cat('\t-NA\t\t\t\tSTR\t(Optional) The NA string. Default: NA\n',file=stderr())
	cat('\t-scale\t\t\t\tSTR\t(Optional) Whether value should be centered, can be "row" or "column" or "none". Default: none\n',file=stderr())
	cat('\t-cluster_rows\t\t\tSTR\t(Optional) If rows should be clustered, can be "F" or "T".  Default: T\n',file=stderr())
	cat('\t-cluster_cols\t\t\tSTR\t(Optional) If columns should be clustered, can be "F" or "T"  Default: T\n',file=stderr())
	cat('\t-clustering_distance_rows\tSTR\t(Optional) Distance measure used in clustering rows. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".  Default: correlation\n',file=stderr())
	cat('\t-clustering_distance_cols\tSTR\t(Optional) Distance measure used in clustering columns. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".  Default: correlation\n',file=stderr())
	cat('\t-clustering_method\t\tSTR\t(Optional) Clustering method used. This must be one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). Default: complete\n',file=stderr())
	cat('\t-treeheight_row\t\t\tINT\t(Optional) The height of a tree for row. Default: 50 points\n',file=stderr())
	cat('\t-treeheight_col\t\t\tINT\t(Optional) The height of a tree for column. Default: 50 points\n',file=stderr())
	cat('\t-show_rownames\t\t\tSTR\t(Optional) If row names are be shown. Default: T\n',file=stderr())
	cat('\t-show_colnames\t\t\tSTR\t(Optional) If column names are be shown. Default: T\n',file=stderr())
	cat('\t-display_numbers\t\tSTR\t(Optional) If the numeric values are also printed to the cells. Default: F\n',file=stderr())
	cat('\t-width\t\t\t\tINT\t(Optional) The width in inches of output. Default: 1024\n',file=stderr())
	cat('\t-height\t\t\t\tINT\t(Optional) The height in inches of output. Default: 768\n',file=stderr())
    cat('\t-color\t\t\t\tSTR\t(Optional) Vector of colors used in heatmap. Default: "blue,darkgreen,gold,red"\n',file=stderr())
	cat('\t-h\t\t\t\t\tShow help\n',file=stderr())
	q(save='no')
}

myPdf='pheatmap.pdf'
myNa='NA'
myMain='Main Title'
myScale='none'
myClusterRows='T'
myClusterCols='T'
myClusterDisRows='correlation'
myClusterDisCols='correlation'
myClusterDisMethod='complete'
myTreeheightRow=50
myTreeheightCol=50
myShowRownames='T'
myShowColnames='T'
myDispalyNumbers='F'
myWidth=1024
myHeight=768
myColor="blue,darkgreen,gold,red"
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
    if(grepl('^-scale?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -scale')
      }else{
        myScale=arg.split[2]
      }
    }
    if(grepl('^-cluster_rows?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -cluster_rows')
      }else{
        myClusterRows=arg.split[2]
      }
    }
    if(grepl('^-cluster_cols?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -cluster_cols')
      }else{
        myClusterCols=arg.split[2]
      }
    }
    if(grepl('^-clustering_distance_rows?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -clustering_distance_rows')
      }else{
        myClusterDisRows=arg.split[2]
      }
    }
    if(grepl('^-clustering_distance_cols?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -clustering_distance_cols')
      }else{
        myClusterDisCols=arg.split[2]
      }
    }
    if(grepl('^-clustering_method?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -clustering_method')
      }else{
        myClusterDisMethod=arg.split[2]
      }
    }
    if(grepl('^-treeheight_row?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
	stop('Please specify the value of -treeheight_row')
      }else{
	myTreeheightRow=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-treeheight_col?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
	stop('Please specify the value of -treeheight_col')
      }else{
	myTreeheightCol=as.numeric(arg.split[2])
      }
    }
    if(grepl('^-show_rownames?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -show_rownames')
      }else{
        myShowRownames=arg.split[2]
      }
    }
    if(grepl('^-show_colnames?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -show_colnames')
      }else{
        myShowColnames=arg.split[2]
      }
    }
    if(grepl('^-display_numbers?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -display_numbers')
      }else{
        myDispalyNumbers=arg.split[2]
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
    if(grepl('^-NA?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -NA')
      }else{
        myNa=arg.split[2]
      }
    }
    if(grepl('^-color?=', arg)){
      arg.split = strsplit(arg, '=', fixed = T)[[1]]
      if(is.na(arg.split[2])){
        stop('Please specify the value of -color')
      }else{
        myColor=arg.split[2]
      }
    }
  }
}

library("pheatmap")
data=read.table(file=myFile,header=TRUE, sep="\t",na.strings=myNa, row.names=1)
dataset=colnames(data)
myLevel1=sapply(strsplit(dataset,"_"),function(a){return(a[length(a)-1])})
myLevel2=sapply(strsplit(dataset,"_"),function(a){return(a[length(a)])})
annotation=data.frame(Var1=myLevel1,Var2=myLevel2)
rownames(annotation)=dataset
typeLevel1=length(levels(annotation$Var1))
typeLevel2=length(levels(annotation$Var2))
Var1=c(1:typeLevel1)
names(Var1)=levels(annotation$Var1)
Var2=c((typeLevel1+1):(typeLevel1+typeLevel2))
names(Var2)=levels(annotation$Var2)
ann_colors = list(Var1 = Var1, Var2=Var2)
#jet.colors <- colorRampPalette(c("blue","darkgreen","gold","red"),bias=1)(100)
myColor=strsplit(myColor,",",fixed=T)
jet.colors <- colorRampPalette(myColor[[1]])(100)
pdf(myPdf)
pheatmap(data,color=jet.colors, annotation = annotation, annotation_colors = ann_colors, scale=myScale, cluster_rows=myClusterRows, cluster_cols=myClusterCols, clustering_distance_rows=myClusterDisRows, clustering_distance_cols=myClusterDisCols,
	clustering_method=myClusterDisMethod, treeheight_row=myTreeheightRow, treeheight_col=myTreeheightCol, show_rownames=as.logical(myShowRownames), show_colnames=as.logical(myShowColnames), display_numbers=myDispalyNumbers, width=myWidth, height=myHeight, main=myMain)

