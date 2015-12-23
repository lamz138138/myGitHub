library("lattice")
library("RColorBrewer")

data=read.table("autosome.phased.dataForPCA",header=T,row.names=1)
dataMatrix=as.matrix(data)

pca=prcomp( t(dataMatrix) )
cmp=rownames(pca$x)
co=sapply(strsplit(cmp,"_"),"[[",2)
or=sapply(strsplit(cmp,"_"),"[[",3)

pcaData=data.frame(
PC1=pca$x[,"PC1"],PC2=pca$x[,"PC2"],
  which="PCA",
  country=co,
  original=or
)

pdf("pca.pdf")

xyplot(
  PC2~PC1 | which, pcaData, cex=2.3,
  pch=c(1:2)[pcaData$country],
  col=c("red","blue","green","black")[pcaData$original],
  key=list(
    space="top",adj=1,
    text=list(levels(pcaData$country),cex=1.3),
    points=list(pch=c(1:2),cex=1.7),
    text=list(levels(pcaData$original),cex=1.3),
    points=list(pch=20,col=c("red","blue","green","black"),cex=1.5),
    rep=FALSE
  )
)

