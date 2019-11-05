#!/bin/env Rscript

library("optparse")

option_list <- list(
  make_option("--outputDir", action="store", type="character", default="output",
              help="Directory of output [default: %default]."),
  make_option("--qValue_enrichment", action="store", type="double", default=0.1,
              help="Thresold of enrichment analysis in kegg [default: %default].")
)

parser <- OptionParser(usage="%prog [options] input", option_list=option_list,
                       description="This script was used to perform enrichment analysis (input are Ensembl IDs).
                       
Following library should be installed: optparse, clusterProfiler, org.Hs.eg.db.")

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
input_file <- args$args[1]

library("clusterProfiler")
library("org.Hs.eg.db")

data=read.delim(file=input_file, header=FALSE)

  
data=data[,1]
data=sapply(data, function(x) sub("\\..*", "", x))

entrezID2symbolEnsembl=bitr(data, fromType = "ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
entrezID2symbolEnsembl=entrezID2symbolEnsembl[!duplicated(entrezID2symbolEnsembl$ENTREZID),]

dir.create(opt$outputDir, showWarnings=FALSE)

getSymbolEnsembl<-function(x){
  y="NA"
  for(i in 1:nrow(x)){
    temp_entrezID=strsplit(x[i]$geneID, "/", fixed=TRUE)
    symbol="NA"
    ensembl="NA"
    for(j in 1:length(temp_entrezID[[1]])){
      symbol_ensembl=entrezID2symbolEnsembl[entrezID2symbolEnsembl$ENTREZID==temp_entrezID[[1]][j],]
      if(j==1){
        symbol=symbol_ensembl$SYMBOL
        ensembl=symbol_ensembl$ENSEMBL
      }else{
        symbol=paste(symbol, symbol_ensembl$SYMBOL, sep="/")
        ensembl=paste(ensembl, symbol_ensembl$ENSEMBL, sep="/")
      }
    }
    temp_output=cbind(x[i], symbol, ensembl)
    temp_output=as.data.frame(temp_output)
    if(i==1){
      y=temp_output
    }else{
      y=rbind(y, temp_output)
    }
  }
  return(y)
} 

# 1. GO
# 1.1. GO_MF
enrich_data=enrichGO(entrezID2symbolEnsembl$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="MF", pAdjustMethod="fdr", qvalueCutoff = opt$qValue_enrichment)
if(!is.null(enrich_data)){
  enrich_data=simplify(enrich_data, cutoff=0.7, by="p.adjust", select_fun = min)
  enrich_final=getSymbolEnsembl(enrich_data)
  enrich_final=enrich_final[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count",
                               "geneID", "symbol", "ensembl")]
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_MF.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_MF.tsv"), sep="\n", append = TRUE)
  write.table(enrich_final, file=file.path(opt$outputDir, "GO_MF.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = TRUE)
}else{
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_MF.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_MF.tsv"), sep="\n", append = TRUE)
}
jpeg(file=file.path(opt$outputDir, "GO_MF.jpg"))
print(dotplot(enrich_data, x="GeneRatio", title="GO_MF"))
dev.off()

# 1.2. GO_BP
enrich_data=enrichGO(entrezID2symbolEnsembl$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP", pAdjustMethod="fdr", qvalueCutoff = opt$qValue_enrichment)
if(!is.null(enrich_data)){
  enrich_data=simplify(enrich_data, cutoff=0.7, by="p.adjust", select_fun = min)
  enrich_final=getSymbolEnsembl(enrich_data)
  enrich_final=enrich_final[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count",
                               "geneID", "symbol", "ensembl")]
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_BP.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_BP.tsv"), sep="\n", append = TRUE)
  write.table(enrich_final, file=file.path(opt$outputDir, "GO_BP.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = TRUE)
}else{
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_BP.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_BP.tsv"), sep="\n", append = TRUE)
}
jpeg(file=file.path(opt$outputDir, "GO_BP.jpg"))
print(dotplot(enrich_data, x="GeneRatio", title="GO_BP"))
dev.off()
  
# 1.3. GO_CC
enrich_data=enrichGO(entrezID2symbolEnsembl$ENTREZID, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC", pAdjustMethod="fdr", qvalueCutoff = opt$qValue_enrichment)
if(!is.null(enrich_data)){
  enrich_data=simplify(enrich_data, cutoff=0.7, by="p.adjust", select_fun = min)
  enrich_final=getSymbolEnsembl(enrich_data)
  enrich_final=enrich_final[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count",
                               "geneID", "symbol", "ensembl")]
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_CC.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_CC.tsv"), sep="\n", append = TRUE)
  write.table(enrich_final, file=file.path(opt$outputDir, "GO_CC.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = TRUE)
}else{
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "GO_CC.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "GO_CC.tsv"), sep="\n", append = TRUE)
}
jpeg(file=file.path(opt$outputDir, "GO_BP.jpg"))
print(dotplot(enrich_data, x="GeneRatio", title="GO_BP"))


# 2. KEGG
enrich_data=enrichKEGG(entrezID2symbolEnsembl$ENTREZID, organism = "hsa", keyType = "ncbi-geneid",  pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = opt$qValue_enrichment, use_internal_data = TRUE)
if(!is.null(enrich_data)){
  enrich_final=getSymbolEnsembl(enrich_data)
  enrich_final=enrich_final[,c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count",
                               "geneID", "symbol", "ensembl")]
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "KEGG.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "KEGG.tsv"), sep="\n", append = TRUE)
  write.table(enrich_final, file=file.path(opt$outputDir, "KEGG.tsv"), sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = TRUE)
}else{
  cat("GeneRatio (k/n): k is the number of genes within that list n, which are annotated to the node, while n is the size of the list of genes of interest.", file=file.path(opt$outputDir, "KEGG.tsv"), sep="\n")
  cat("BgRatio (M/N): M is  the number of genes within that distribution that are annotated (either directly or indirectly) to the node of interest, while N is the total number of genes in the background distribution (universe).", file=file.path(opt$outputDir, "KEGG.tsv"), sep="\n", append = TRUE)
}
jpeg(file=file.path(opt$outputDir, "KEGG.jpg"))
print(dotplot(enrich_data, x="GeneRatio", title="KEGG"))
dev.off()

# 3. DAVID
#enrich_data=enrichDAVID(entrezID2symbolEnsembl$ENTREZID, idType = "ENTREZ_GENE_ID", pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = opt$qValue_enrichment)

