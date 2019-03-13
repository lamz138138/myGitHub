#!/bin/bash

usage(){
 cat << EOF
Description:
    This script was used to get different expressed genes with cuffdiff or DESeq2.
    Please use -t to indicate the directory contain different-expression-DESeq2.R
Usage:
    different-expression.sh 
Options:
    -o Output directory [default: output]
    -t Dirctory contain different-expression-DESeq2.R [default: PWD]
    -c Use cuffdiff to find different expressed genes
    -d Use DESeq2.R to find different expressed genes
    -r Reference in gtf format
    -t Number of threads [default: 1]
    -l Library type (one of fr-unstranded, fr-firststrand, fr-secondstrand) [default: fr-unstranded]
EOF
    exit 0
}

[ $# -eq 0 ] && usage
[ -z $refPath ] && echo "Please assign -r" && exit 1

outputDir=$PWD/output
toolDir=$PWD
threadNum=1
libraryType="fr-unstranded"
while getopts "ho:t:cdr:" OPTION
 do
  case $OPTION in
        h) usage;;
        o) outputDir=$OPTARG;;
        t) toolDir=$OPTARG;;
        c) useCuffdiff=1;;
        d) useDESeq2=1;;
        r) refPath=$OPTARG;;
        t) threadNum=$OPTARG;;
        l) libraryType=$OPTARG;;
  esac
 done
shift $((OPTIND - 1))
 

outputDir=$(readlink -f $outputDir)
[ -e $outputDir ] && mv $outputDir ${outputDir}_bak
mkdir -p $outputDir/log
cd $outputDir 
 
if [ ! -z $useCuffdiff ]; then
	cuffcompare $refPath $refPath
	cuffdiff 	-p $threadNum -no-update-check -emit-count-tables cuffcmp.combined.gtf
else

fi 
 
 
# 1. get co-expression and perform enrichment analysis for first gene ($1 is ensembl ID of gene, $3 is file containning gene expression in tumor samples,
# $4 is file containning gene expression in normal sample, and $5 is file containniing entrezID of transcription fator genes)
# orrelation-heatmap.R --gene ENSG00000156273.14 expression-tumor-brca-protein.tsv --tf tfs-checkpoint-human.tsv  --outputDir output 
$toolDir/correlation-heatmap.R --gene $1 $3 $4 --tf $5 --outputDir $outputDir/output-first >$outputDir/log/correlation-heatmap-first.log 2>$outputDir/log/correlation-heatmap-first.err

# 2. get co-expression and perform enrichment analysis for second gene (tf) ($1 is ensembl ID of gene, $3 is file containning gene expression in tumor samples,
# $4 is file containning gene expression in normal sample, and $5 is file containniing entrezID of transcription fator genes)
# orrelation-heatmap.R --gene ENSG00000156273.14 expression-tumor-brca-protein.tsv --tf tfs-checkpoint-human.tsv  --outputDir output 
$toolDir/correlation-heatmap.R --gene $2 $3 $4 --tf $5 --outputDir $outputDir/output-second >$outputDir/log/correlation-heatmap-second.log 2>$outputDir/log/correlation-heatmap-second.err

# 3. get overlap of coexpress gene between first and second gene
$toolDir/selectionOutput.pl --colmun11 6 --colmun12 6 $outputDir/output-first/correlation-filtered-pearson.tsv $outputDir/output-second/correlation-filtered-pearson.tsv 2>/dev/null | cut -f6 >$outputDir/overlap-pearson.tsv
$toolDir/selectionOutput.pl --colmun11 6 --colmun12 6 $outputDir/output-first/correlation-filtered-spearman.tsv $outputDir/output-second/correlation-filtered-spearman.tsv 2>/dev/null | cut -f6 >$outputDir/overlap-spearman.tsv

# 4. enrichment analysis of overlap genes
$toolDir/enrichment-clusterProfiler.R $outputDir/overlap-pearson.tsv --tf $5  --outputDir $outputDir/overlap-enrichment --cor_method pearson >$outputDir/log/overlap-enrichment-pearson.log 2>$outputDir/log/overlap-enrichment-pearson.err
$toolDir/enrichment-clusterProfiler.R $outputDir/overlap-spearman.tsv --tf $5  --outputDir $outputDir/overlap-enrichment --cor_method spearman >$outputDir/log/overlap-enrichment-spearman.log 2>$outputDir/log/overlap-enrichment-spearman.err

