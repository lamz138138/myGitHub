#!/bin/env Rscript

# please refer to http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html for algorithm

library("optparse")

option_list <- list(
  make_option("--snpInRegion", action="store", type="integer", default=1,
              help="snps in candidate region [default: %default]."),
  make_option("--totalSNPs", action="store", type="integer", default=1,
              help="total snps [default: %default]."),
  make_option("--regionSize", action="store", type="integer", default=1,
              help="size of candidate region [default: %default]."),
  make_option("--genomeSize", action="store", type="integer", default=1,
              help="size of genome [default: %default].")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list,
                       description="This script was used to test whether candidate enrich snps by using phyper.")

args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

q=opt$snpInRegion  # snps in candidate region
m=opt$totalSNPs # total snps 
n=opt$genomeSize - opt$totalSNPs # total non-snps or genome size - total snps
k=opt$regionSize # candidate region size

phyper(q-1, m, n, k, lower.tail=FALSE)

