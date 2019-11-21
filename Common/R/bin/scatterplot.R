#!/bin/env Rscript

library("optparse")


option_list <- list(
  make_option(c("-o", "--outputDir"), action="store", type="character", default="output",
              help="Directory of output [default: %default]."),
  make_option(c("-n", "--outputName"), action="store", type="character", default="scatter",
              help="Directory of output [default: %default]."),
  make_option(c("-x", "--xColumn"), action="store", type="character", default="output",
              help="Column to be used for X"),
  make_option(c("-y", "--yColumn"), action="store", type="character", default="output",
              help="Column to be used for y."),
  make_option(c("-l", "--log"), action="store_true", default=FALSE,
              help="Log transformation (xTrans and yTrans must be given)"),
  make_option(c("--xTrans"), action="store", type="character", default="log10",
              help="Log transformation for X (--log must be given, log2, log10, sqrt.....) [default: %default]."),
  make_option(c("--yTrans"), action="store", type="character", default="log10",
              help="Log transformation for Y (--log must be given, log2, log10, sqrt.....) [default: %default]."),
  make_option(c("--xIntercept"), action="store", type="character", default="NA",
              help="Parameters that control the position of the line in X [default: %default]."),
  make_option(c("--yIntercept"), action="store", type="character", default="NA",
              help="Parameters that control the position of the line in Y [default: %default]."),
  make_option(c("--abline"), action="store_true", default=FALSE,
              help="Draw X=Y")
)

parser <- OptionParser(usage="%prog [options] input", option_list=option_list,
                       description="This script was used to draw scatter plot (must contain header).")

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options

library("ggplot2")
input <- args$args[1]

data=read.delim(file=input,row.names = 1, header=TRUE, check.names = FALSE)

g<-ggplot(data, aes_string(x=opt$xColumn, y=opt$yColumn)) + geom_point() + expand_limits(x = 0, y = 0)

if(opt$log){
  g<-g+scale_x_continuous(trans=opt$xTrans)+scale_y_continuous(trans=opt$yTrans)
}

if(opt$xIntercept!="NA"){
  g<-g+geom_vline(xintercept = as.double(opt$xIntercept), color="red")
}

if(opt$yIntercept!="NA"){
  g<-g+geom_hline(yintercept = as.double(opt$yIntercept), color="red")
}
  
if(opt$abline){
  g<-g+geom_abline(color="red")
}
  
jpeg(file=file.path(opt$outputDir, paste(opt$outputName, ".jpg", sep="")))
print(g)
#dev.off()

# color point and solve overlap of point: geom_point(aes(colour=Type), position = position_jitter(h=0.05,w=0.15))
# text according to criterial (Normal, Tumor and Symbol are column name): geom_text(aes(label=ifelse(Normal==0 & Tumor>2,as.character(Symbol),'')), position = position_jitter(h=0.05,w=0.15), hjust=0, vjust=0,size=2) 
