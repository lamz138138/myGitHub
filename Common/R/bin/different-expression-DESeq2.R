#!/bin/env Rscript

# 2019-09-20
# developed by zhongxm

# countData
# SC5-LA  SC5-LV  SC5-RA  SC5-RV  SC6-LA  SC6-LV  SC6-RA  SC6-RV  SC8-LA  SC8-LV  SC8-RA  SC8-RV
#A1BG    107     59      109     35      149     62      180     32      139     88      106     104
#A1BG-AS1        16      19      24      19      68      47      111     50      42      30      35      24
#A1CF    77      66      78      85      13      13      10      14      12      17      16      10
#A2M     534     611     663     547     534     416     553     370     388     494     386     258
#A2M-AS1 5475    4467    4152    4949    8164    5142    5204    4629    8150    10210   5940    10030
#A2ML1   22      10      18      24      25      20      26      20      17      38      23      8

# metaData
#sample	tissue	individual
#SC5-LA	A	SC5
#SC5-LV	V	SC5
#SC5-RA	A	SC5
#SC5-RV	V	SC5
#SC6-LA	A	SC6
#SC6-LV	V	SC6
#SC6-RA	A	SC6
#SC6-RV	V	SC6
#SC8-LA	A	SC8
#SC8-LV	V	SC8
#SC8-RA	A	SC8
#SC8-RV	V	SC8

library("optparse")
library("DESeq2")



option_list <- list(
  make_option("--outputDir", action="store", type="character", default="output",
              help="Directory of output [default: %default].")
)

parser <- OptionParser(usage="%prog [options] input_count input_meta", option_list=option_list,
                       description="This script was used to get different expresse genes.")

args <- parse_args(parser, positional_arguments = 2)
opt <- args$options
input_count <- args$args[1]
input_meta <- args$args[2]


countData=read.csv(file=input_count, header=TRUE, sep="\t")
metaData=read.csv(file=input_meta, header=TRUE, sep="\t") # header is sample, tissue, individual

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~tissue, tidy = TRUE)

dds <- DESeq(dds)

#res <- results(dds, contrast=c("tissue", "LA", "LV"))
write.table(results(dds, contrast=c("tissue", "A", "V"), tidy=TRUE), file=file.path(opt$outputDir, "result-AV.tsv"), sep="\t", quote=FALSE, row.names=FALSE)