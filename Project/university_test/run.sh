#!/bin/bash

gunzip test.assignment.txt.gz
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/{print "chr"$1,$2-1,$2;}' test.assignment.txt >snp.bed3
sort -k1,1 -k2,2n snp.bed3 >snp_sorted.bed3
rm snp.bed3 
echo "total snps"
wc -l snp_sorted.bed3 | cut -f1 -d " " # 970572
echo "genome size (only chr1~chr22,chrX,chrY,chrMT)"
awk 'BEGIN{FS=OFS="\t"} $1!~/_/{sum+=$2;} END{print sum;}' /mnt/share/share/data/chr.size/hg19.size # 3095693983


# 1. repeat_makser
mkdir repeat_masker && cd repeat_masker
# 1) get chromosome region of LINE
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql
zcat rmsk.txt.gz | awk 'BEGIN{FS=OFS="\t"} $6!~/_/ && $12=="LINE"{print $6,$7,$8,"name","0",$10;}' >LINE.bed6
sort -k1,1 -k2,2n LINE.bed6 >LINE_sorted.bed6
bedtools merge -i LINE_sorted.bed6 >LINE_merge.bed3
rm LINE.bed6
# 2) get overlap between snp and LINE
sort -k1,1 -k2,2n LINE_merge.bed3 >LINE_merge_sorted.bed3
rm LINE_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b LINE_merge_sorted.bed3 >intersect.tsv
echo "total snps in LINE region"
wc -l intersect.tsv | cut -f1 -d " " # 175485
echo "LINE size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' LINE_merge_sorted.bed3 # 630587623
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed LINE_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} {GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.3690
# 4) enrichment test
./phyper.R --snpInRegion 175485 --totalSNPs 970572 --regionSize 630587623 --genomeSize 3095693983


# 2. HAIB_Methyl
mkdir HAIB_Methyl && cd HAIB_Methyl
# 1) get chromosome region of methylation
wget -r -np -nd --accept=bed.gz -e robots=off http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/
zcat *bed.gz | awk 'BEGIN{FS=OFS="\t"} $0!~/^track/{print $1,$2,$3,$4,$5,$6;}' > methy.bed6
sort -k1,1 -k2,2n methy.bed6 >methy_sorted.bed6
bedtools merge -i methy_sorted.bed6 >methy_merge.bed3
rm methy.bed6
# 2) get overlap between snp and methy
sort -k1,1 -k2,2n methy_merge.bed3 >methy_merge_sorted.bed3
rm methy_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b methy_merge_sorted.bed3 >intersect.tsv
echo "total snps in methy region"
wc -l intersect.tsv | cut -f1 -d " " # 4876
echo "methy size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' methy_merge_sorted.bed3 # 2724380
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed methy_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} {GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 1.0000
# 4) enrichment test
./phyper.R --snpInRegion 4876 --totalSNPs 970572 --regionSize 2724380 --genomeSize 3095693983


# 3. enhancer
mkdir enhancer && cd enhancer
# 1) get chromosome region of enhancer
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vistaEnhancers.txt.gz
zcat *txt.gz | cut -f2-6 > enhancer.bed5
sort -k1,1 -k2,2n enhancer.bed5 >enhancer_sorted.bed5
bedtools merge -i enhancer_sorted.bed5 >enhancer_merge.bed3
rm enhancer.bed5
# 2) get overlap between snp and enhancer
sort -k1,1 -k2,2n enhancer_merge.bed3 >enhancer_merge_sorted.bed3
rm enhancer_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b enhancer_merge_sorted.bed3 >intersect.tsv
echo "total snps in enhancer region"
wc -l intersect.tsv | cut -f1 -d " " # 683
echo "enhancer size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' enhancer_merge_sorted.bed3 # 2027471
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed enhancer_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} {GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.4003
# 4) enrichment test
./phyper.R --snpInRegion 683 --totalSNPs 970572 --regionSize 2027471 --genomeSize 3095693983


# 4. TSS
mkdir TSS && cd TSS
# 1) get chromosome region of TSS
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/switchDbTss.txt.gz
zcat *txt.gz | cut -f2-7 > tss.bed6
sort -k1,1 -k2,2n tss.bed6 >tss_sorted.bed6
bedtools merge -i tss_sorted.bed6 >tss_merge.bed3
rm tss.bed6
# 2) get overlap between snp and tss
sort -k1,1 -k2,2n tss_merge.bed3 >tss_merge_sorted.bed3
rm tss_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b tss_merge_sorted.bed3 >intersect.tsv
echo "total snps in tss region"
wc -l intersect.tsv | cut -f1 -d " " # 48
echo "tss size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' tss_merge_sorted.bed3 # 131748
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed tss_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} {GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.5075
# 4) enrichment test
./phyper.R --snpInRegion 48 --totalSNPs 970572 --regionSize 131748 --genomeSize 3095693983


# 5. gene
mkdir gene && cd gene
# 1) get chromosome region of cancer census gene
sftp "zhongxm@pku.edu.cn"@sftp-cancer.sanger.ac.uk
get /files/grch37/cosmic/v84/cancer_gene_census.csv
cp cancer_gene_census.csv cancer_gene_census.csv.bak
#sed -i "s/.*,\(.*:[0-9]+-[0-9]+\),*/\1/" cancer_gene_census.csv
sed -i "s/.*,\([0-9XY]*:[0-9]*-[0-9]*\),.*/\1/" cancer_gene_census.csv
# manually revised gene.bed3 for duplication gene
sed -i "/1585451/d;/1331527/d" cancer_gene_census.csv
awk 'BEGIN{FS=","; OFS="\t"} NR>1{chr=$1;sub(":.*","",chr);chr="chr"chr; start=$1; sub("-.*","",start); sub(".*:","",start); end=$1; sub(".*-","",end); print chr,start,end;}' cancer_gene_census.csv >gene.bed3
echo -e "chrX\t1314894\t1331527" >>gene.bed3
echo -e "chrY\t1314894\t1331527" >>gene.bed3
echo -e "chrX\t1584372\t1585451" >>gene.bed3
echo -e "chrY\t1584372\t1585451" >>gene.bed3
sort -k1,1 -k2,2n gene.bed3 >gene_sorted.bed3
awk 'BEGIN{FS=OFS="\t"} $2!=""' gene_sorted.bed3 >gene_sorted.bed3.v2
mv gene_sorted.bed3.v2 gene_sorted.bed3
bedtools merge -i gene_sorted.bed3 >gene_merge.bed3
rm gene.bed3
# 2) get overlap between snp and gene
sort -k1,1 -k2,2n gene_merge.bed3 >gene_merge_sorted.bed3
rm gene_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b gene_merge_sorted.bed3 >intersect.tsv
echo "total snps in gene region"
wc -l intersect.tsv | cut -f1 -d " " # 23524
echo "gene size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' gene_merge_sorted.bed3 # 67627758
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed gene_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} {GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.4096
# 4) enrichment test
./phyper.R --snpInRegion 23524 --totalSNPs 970572 --regionSize 67627758 --genomeSize 3095693983

