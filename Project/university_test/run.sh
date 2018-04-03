#!/bin/bash

# snp == mutation in this script
# bed2gpe.pl, bedFeature.pl and gpeFeature.pl is inhouse script

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
<<mark
# 1) get chromosome region of simple_repeat
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql
zcat rmsk.txt.gz | awk 'BEGIN{FS=OFS="\t"} $6!~/_/ && $12=="Simple_repeat"{print $6,$7,$8,"name","0",$10;}' >simple_repeat.bed6
sort -k1,1 -k2,2n simple_repeat.bed6 >simple_repeat_sorted.bed6
bedtools merge -i simple_repeat_sorted.bed6 >simple_repeat_merge.bed3
rm simple_repeat.bed6
# 2) get overlap between snp and simple_repeat
sort -k1,1 -k2,2n simple_repeat_merge.bed3 >simple_repeat_merge_sorted.bed3
rm simple_repeat_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b simple_repeat_merge_sorted.bed3 >intersect.tsv
echo "total snps in simple_repeat region"
wc -l intersect.tsv | cut -f1 -d " " # 70313
echo "simple_repeat size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' simple_repeat_merge_sorted.bed3 # 25984317
# 3) GC content
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed simple_repeat_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.3212
# 4) enrichment test
./phyper.R --snpInRegion 70313 --totalSNPs 970572 --regionSize 25984317 --genomeSize 3095693983
mark

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
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed LINE_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.3690
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
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed methy_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 1.0000
# 4) enrichment test
./phyper.R --snpInRegion 4876 --totalSNPs 970572 --regionSize 2724380 --genomeSize 3095693983
# 5) mutation type of enriched snp
awk 'BEGIN{FS=OFS="\t"} {print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >snp.bed3+
sort -k1,1 -k2,2n snp.bed3+ >snp_sorted.bed3+
rm snp.bed3+
bedtools intersect -wa -wb -a intersect.tsv -b snp_sorted.bed3+ | cut -f1-3,10-12 >snp_type.bed3+
echo "mutation type and count"
cut -f4,5 snp_type.bed3+ | sort | uniq -c
# whether downstream base of C is G
awk 'BEGIN{FS=OFS="\t"} $4=="C" && $5=="T"{print $1,$2+1,$3+1;}' snp_type.bed3+ >downstream_C.bed3
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed downstream_C.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{A+=$6; C+=$7; G+=$8; T+=$9; sum+=$12;} END{print A; print C; print G; print T; print sum; printf("%.4f", A/sum); \
	printf("%.4f", C/sum); printf("%.4f", G/sum); printf("%.4f", T/sum);}' # 1.00
# whether upstream base of G is C
awk 'BEGIN{FS=OFS="\t"} $4=="G" && $5=="A"{print $1,$2-1,$3-1;}' snp_type.bed3+ >upstream_G.bed3
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed upstream_G.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{A+=$6; C+=$7; G+=$8; T+=$9; sum+=$12;} END{{print A; print C; print G; print T; print sum; printf("%.4f", A/sum); \
	printf("%.4f", C/sum); printf("%.4f", G/sum); printf("%.4f", T/sum);}' # 1.00 


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
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed enhancer_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.4003
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
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed tss_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.5074
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
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed gene_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.4096
# 4) enrichment test
./phyper.R --snpInRegion 23524 --totalSNPs 970572 --regionSize 67627758 --genomeSize 3095693983
# 5) enrich in GO, DAVID 
cp cancer_gene_census.csv.bak cancer_gene_census.csv.v2
sed -i "s/.*,\([0-9]*,[0-9XY]*:[0-9]*-[0-9]*\),.*/\1/" cancer_gene_census.csv.v2
sed -i "/1585451/d;/1331527/d" cancer_gene_census.csv.v2
awk 'BEGIN{FS=","; OFS="\t"} NR>1{chr=$2;sub(":.*","",chr);chr="chr"chr; start=$2; sub("-.*","",start); sub(".*:","",start); end=$2; sub(".*-","",end); print chr,start,end,$1;}' cancer_gene_census.csv.v2 >gene_entrezID.bed3+
echo -e "chrX\t1314894\t1331527\t64109" >>gene_entrezID.bed3+
echo -e "chrY\t1314894\t1331527\t64109" >>gene_entrezID.bed3+
echo -e "chrX\t1584372\t1585451\t286530" >>gene_entrezID.bed3+
echo -e "chrY\t1584372\t1585451\t286530" >>gene_entrezID.bed3+
rm cancer_gene_census.csv.v2
sort -k1,1 -k2,2n gene_entrezID.bed3+ >gene_entrezID_sorted.bed3+
awk 'BEGIN{FS=OFS="\t"} $2!=""' gene_entrezID_sorted.bed3+ >gene_entrezID_sorted.bed3+.v2
mv gene_entrezID_sorted.bed3+.v2 gene_entrezID_sorted.bed3+
rm gene_entrezID.bed3+ 
bedtools intersect -wa -wb -a intersect.tsv -b gene_entrezID_sorted.bed3+  | cut -f1-3,7-10 | sort | uniq >entrezID.bed3+
# 6) exon, intron, cds, utr, intergenic
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
zcat gene2refseq.gz | awk 'BEGIN{FS=OFS="\t"} $1=="9606" && $4!="-"{gsub("\\..*","",$4); print $2,$4;}' | sort | uniq > entrezID_refSeq.tsv
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz
selectionOutput.pl entrezID.bed3+ entrezID_refSeq.tsv >entrezID_refSeq_gene.tsv 2>/dev/null
selectionOutput.pl --colmun11 2 --colmun12 2 entrezID_refSeq_gene.tsv refGene.txt >refGene_select.gpe 2>/dev/null
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i refGene_select.gpe >intron.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e refGene_select.gpe >exon.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -u refGene_select.gpe >utr.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -c refGene_select.gpe >cds.bed
gpeFeature.pl -b --upstream 1000 --downstream 1000 -g /mnt/share/share/data/chr.size/hg19.size refGene_select.gpe >intergenic.bed
gpe2bed.pl -b refGene_select.gpe >refGene_select.bed
bedFeature.pl -s refGene_select.bed >splice_site.bed
for type in intron exon utr cds intergenic splice_site
  do
    sort -k1,1 -k2,2n "$type".bed >"$type"_sorted.bed
    rm "$type".bed
    bedtools merge -i "$type"_sorted.bed >"$type"_merge.bed3
    sort -k1,1 -k2,2n "$type"_merge.bed3 >"$type"_merge_sorted.bed3
    rm "$type"_merge.bed3
    bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b "$type"_merge_sorted.bed3 >intersect_"$type".tsv
    echo "total snps in $type region"
    wc -l intersect_"$type".tsv | cut -f1 -d " " 
    echo "$type size"
    awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' "$type"_merge_sorted.bed3
  done
# type (snp, size): intron (28164, 81166462), exon (1377, 3849024), utr (630, 2112623), cds (827, 1955596), intergenic (576, 1821835), splice_site (0, 3900)
for type in intron exon utr cds intergenic splice_site
do 
echo "$type GC"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed "$type"_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' 
done
# intron (0.4076), exon (0.4776), utr (0.4482), cds (0.5127), intergenic (0.5032), splice_site (0.5226)


# 6. GC content
mkdir gc && cd gc
# 1) v1
awk 'BEGIN{FS=OFS="\t"} $1!~/_/{myIndex=int($2/1000000); print $1,"0","1000000"; for(i=1; i<myIndex; i++){print $1, i*1000000+1, (i+1)*1000000;} }' /mnt/share/share/data/chr.size/hg19.size >region.bed3
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed region.bed3 | cut -f1-3,5 >region_gc.bed3 # chrMT is skiped
<<mark
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4<=0.25' region_gc.bed3 >region_1.bed3+ # GC in [0,0.25] -> 231
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.25 && $4<=0.5' region_gc.bed3 >region_2.bed3+ # GC in (0.25,0.50] -> 2753
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.5 && $4<=0.75' region_gc.bed3 >region_3.bed3+ # GC in (0.50,0.75] -> 105
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.75 && $4<=0.1' region_gc.bed3 >region_4.bed3+ # GC in (0.75,1.00] -> 0
mark
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.35 && $4<=0.4' region_gc.bed3 >region_1.bed3+ # GC in (0.35,0.4] -> 1281
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.4 && $4<=0.45' region_gc.bed3 >region_1.bed3+ # GC in (0.4,0.45] -> 955
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>0.45 && $4<=0.5' region_gc.bed3 >region_1.bed3+ # GC in (0.45,0.5] -> 355
sort -k1,1 -k2,2n region_1.bed3+ >region_1_sorted.bed3+
sort -k1,1 -k2,2n region_2.bed3+ >region_2_sorted.bed3+
sort -k1,1 -k2,2n region_3.bed3+ >region_3_sorted.bed3+
rm region_1.bed3+ region_2.bed3+ region_3.bed3+
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b region_1_sorted.bed3+ >intersect_1.tsv
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b region_2_sorted.bed3+ >intersect_2.tsv
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b region_2_sorted.bed3+ >intersect_2.tsv
cut -f4-6 intersect_1.tsv | sort | uniq -c | sed "s/\s*//" | cut -d " " -f1 >snp_count_1.tsv
cut -f4-6 intersect_2.tsv | sort | uniq -c | sed "s/\s*//" | cut -d " " -f1 >snp_count_2.tsv
cut -f4-6 intersect_3.tsv | sort | uniq -c | sed "s/\s*//" | cut -d " " -f1 >snp_count_3.tsv 
wilcox.R -file1=snp_count_1.tsv -file2=snp_count_2.tsv # 0.0003109803
wilcox.R -file1=snp_count_1.tsv -file2=snp_count_3.tsv # 1.199875e-06
wilcox.R -file1=snp_count_2.tsv -file2=snp_count_3.tsv # 3.763976e-18
# 2) v2
perl /rd1/user/zhongxm/bin/selectionOutput.pl --colmun12 2 <(grep -v "^#" ../snpEff_annotation/snpEff_genes.txt | cut -f3) ../gene/refGene.txt >refGene_longest.gpe 2>/dev/null
awk 'BEGIN{FS=OFS="\t"} {print $6-$5;}' refGene_longest.gpe >refGene_longest_length.txt
sed -i "1i #head" refGene_longest_length.txt
boxplot.R -p=refGene_longest_length.pdf -file=refGene_longest_length.txt
awk 'BEGIN{FS=OFS="\t"} $1!~/_/{myIndex=int($2/10000); print $1,"0","10000"; for(i=1; i<myIndex; i++){print $1, i*10000+1, (i+1)*10000;} }' /mnt/share/share/data/chr.size/hg19.size >region.bed
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed region.bed | cut -f1-3,5 >region_gc.bed # chrMT is skiped
boxplot.R -p=gc.pdf -column=4 -file=region_gc.bed
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>=0 && $4<0.25' region_gc.bed >region_gc_0_0.25.bed3+ # GC in [0,0.25) -> 
awk 'BEGIN{FS=OFS="\t"} NR>1 && $4>=0.55 && $4<=1' region_gc.bed >region_gc_0.55_1.bed3+ # GC in [0.55,1] -> 
for i in 0.25 0.3 0.35 0.4 0.45 0.5
  do
    j=$( echo "$i+0.05" | bc )
    awk -v i=$i 'BEGIN{FS=OFS="\t"} NR>1 && $4>=i && $4<i+0.05' region_gc.bed >region_gc_"$i"_0"$j".bed3+
  done
ls region_gc_0* | while read file
  do
    sample=$( echo $file | sed 's/gc/gc_sorted/' )
    sort -k1,1 -k2,2n $file >$sample
    rm $file
    intersect=$( echo $file | sed 's/region_gc/intersect/;s/bed3+/tsv/' )
    bedtools intersect -wa -wb -a ../bmr/snp_sorted.bed3 -b $sample >"$intersect"
    snpCount=$( echo $file | sed 's/region/snp_count/;s/bed3+/tsv/' )
    cut -f4-6 $intersect | sort | uniq -c | sed "s/\s*//" | cut -d " " -f1 >$snpCount
  done
export tempArrary=("snp_count_gc_0_0.25.tsv" "snp_count_gc_0.25_0.30.tsv" "snp_count_gc_0.3_0.35.tsv" "snp_count_gc_0.35_0.40.tsv" "snp_count_gc_0.4_0.45.tsv" "snp_count_gc_0.45_0.50.tsv" "snp_count_gc_0.5_0.55.tsv" "snp_count_gc_0.55_1.tsv")
for((i=0; i<7; i++))
  do
    for((j=$i+1; j<8; j++))
      do
	echo ${tempArrary["$i"]} ${tempArrary["$j"]}
	wilcox.R -file1=${tempArrary["$i"]} -file2=${tempArrary["$j"]}
      done
  done


# 7. exon, intron, intergenic, cds
mkdir exon_intron_intergenic_cds && cd exon_intron_intergenic_cds
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i refGene.txt >intron.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e refGene.txt >exon.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -u refGene.txt >utr.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -c refGene.txt >cds.bed
for type in intron exon utr cds
  do
    sort -k1,1 -k2,2n "$type".bed >"$type"_sorted.bed
    rm "$type".bed
    bedtools merge -i "$type"_sorted.bed >"$type"_merge.bed3
    sort -k1,1 -k2,2n "$type"_merge.bed3 >"$type"_merge_sorted.bed3
    rm "$type"_merge.bed3
    bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b "$type"_merge_sorted.bed3 >intersect_"$type".tsv
    echo "total snps in $type region"
    wc -l intersect_"$type".tsv | cut -f1 -d " " 
    echo "$type size"
    awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' "$type"_merge_sorted.bed3
  done
cat *merge_sorted.bed3 >non_inter.bed3
sort -k1,1 -k2,2n non_inter.bed3 >non_inter_sorted.bed3
rm non_inter.bed3
bedtools merge -i non_inter_sorted.bed3 >non_merge.bed3
bedtools complement -i non_merge.bed3 -g /mnt/share/share/data/chr.size/hg19.size >intergenic.bed
rm non*
for type in intergenic
  do
    sort -k1,1 -k2,2n "$type".bed >"$type"_sorted.bed
    rm "$type".bed
    bedtools merge -i "$type"_sorted.bed >"$type"_merge.bed3
    sort -k1,1 -k2,2n "$type"_merge.bed3 >"$type"_merge_sorted.bed3
    rm "$type"_merge.bed3
    bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b "$type"_merge_sorted.bed3 >intersect_"$type".tsv
    echo "total snps in $type region"
    wc -l intersect_"$type".tsv | cut -f1 -d " " 
    echo "$type size"
    awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' "$type"_merge_sorted.bed3
  done
# type (snp, size): intron (418119, 1300175136), exon (26580, 84102179), utr (15851, 53601847), cds (12370, 35648926), intergenic (530662, 1769068317)
for type in intron exon utr cds intergenic
do 
echo "$type GC"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed "$type"_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' 
done
# intron (0.4121), exon (0.4868), utr (0.4645), cds (0.5247), intergenic (0.3490)


# 8. splice site
mkdir splice_site && cd splice_site
ln -s ../exon_intron_intergenic_cds/refGene.txt .
gpe2bed.pl -b refGene.txt >refGene.bed12
bedFeature.pl -s refGene.bed12 >splice_site.bed6+
sort -k1,1 -k2,2n splice_site.bed6+ >splice_site_sorted.bed6+
rm splice_site.bed6+
bedtools merge -i splice_site_sorted.bed6+ >splice_site_merge.bed3
sort -k1,1 -k2,2n splice_site_merge.bed3 >splice_site_merge_sorted.bed3
rm splice_site_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b splice_site_merge_sorted.bed3 >intersect.tsv
echo "total snps in splice sites region"
wc -l intersect.tsv | cut -f1 -d " " # 32
echo "splice sites size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' splice_site_merge_sorted.bed3 # 144385 
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed splice_site_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.5079


# 9. miRNA
mkdir miRNA && cd miRNA
# miRNA targe were inhouse data predicted by miRanda and targetScan
cut -f22-24 miRNA_miRanda.tsv miRNA_targetScan.tsv >miRNA_target.bed3
sort -k1,1 -k2,2n miRNA_target.bed3 >miRNA_target_sorted.bed3
awk 'BEGIN{FS=OFS="\t"} $2!=""' miRNA_target_sorted.bed3 >miRNA_target_sorted.bed3.v2
mv miRNA_target_sorted.bed3.v2 miRNA_target_sorted.bed3
rm miRNA_target.bed3
bedtools merge -i miRNA_target_sorted.bed3 >miRNA_target_merge.bed3
sort -k1,1 -k2,2n miRNA_target_merge.bed3 >miRNA_target_merge_sorted.bed3
rm miRNA_target_merge.bed3
bedtools intersect -wa -wb -a ../snp_sorted.bed3 -b miRNA_target_merge_sorted.bed3 >intersect.tsv
echo "total snps in miRNA target region"
wc -l intersect.tsv | cut -f1 -d " " # 34490
echo "miRNA target size"
awk 'BEGIN{FS=OFS="\t"} {sum+=$3-$2;} END{print sum;}' miRNA_target_merge_sorted.bed3 # 115560393
echo "output GC bases, total bases, GC percentage"
bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed miRNA_target_merge_sorted.bed3 | awk 'BEGIN{FS=OFS="\t"} NR>1{GC+=$7+$8; sum+=$12;} END{print GC; print sum; printf("%.4f", GC/sum);}' # 0.4126


# 10. mutation overlap
mkdir mutation_overlap && cd mutation_overlap
awk 'BEGIN{FS=OFS="\t"} $5=="TCGA-AA-3516-01A-02D-1554-10" && $1!~/GL/{print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >sample_1.bed3+
awk 'BEGIN{FS=OFS="\t"} $5=="TCGA-AA-A01R-01A-21D-A17O-10" && $1!~/GL/{print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >sample_2.bed3+
awk 'BEGIN{FS=OFS="\t"} $5=="TCGA-AD-6964-01A-11D-1924-10" && $1!~/GL/{print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >sample_3.bed3+
awk 'BEGIN{FS=OFS="\t"} $5=="TCGA-AD-A5EJ-01A-11D-A28G-10" && $1!~/GL/{print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >sample_4.bed3+
awk 'BEGIN{FS=OFS="\t"} $5=="TCGA-AZ-6601-01A-11D-1771-10" && $1!~/GL/{print "chr"$1,$2-1,$2,$3,$4,$5;}' ../test.assignment.txt >sample_5.bed3+
sort -k1,1 -k2,2n sample_1.bed3+ >sample_1_sorted.bed3+
sort -k1,1 -k2,2n sample_2.bed3+ >sample_2_sorted.bed3+
sort -k1,1 -k2,2n sample_3.bed3+ >sample_3_sorted.bed3+
sort -k1,1 -k2,2n sample_4.bed3+ >sample_4_sorted.bed3+
sort -k1,1 -k2,2n sample_5.bed3+ >sample_5_sorted.bed3+
bedtools intersect -wa -wb -a sample_1_sorted.bed3+ -b sample_2_sorted.bed3+ >intersect_1_2.tsv
bedtools intersect -wa -wb -a sample_1_sorted.bed3+ -b sample_3_sorted.bed3+ >intersect_1_3.tsv
bedtools intersect -wa -wb -a sample_1_sorted.bed3+ -b sample_4_sorted.bed3+ >intersect_1_4.tsv
bedtools intersect -wa -wb -a sample_1_sorted.bed3+ -b sample_5_sorted.bed3+ >intersect_1_5.tsv
bedtools intersect -wa -wb -a sample_2_sorted.bed3+ -b sample_3_sorted.bed3+ >intersect_2_3.tsv
bedtools intersect -wa -wb -a sample_2_sorted.bed3+ -b sample_4_sorted.bed3+ >intersect_2_4.tsv
bedtools intersect -wa -wb -a sample_2_sorted.bed3+ -b sample_5_sorted.bed3+ >intersect_2_5.tsv
bedtools intersect -wa -wb -a sample_3_sorted.bed3+ -b sample_4_sorted.bed3+ >intersect_3_4.tsv
bedtools intersect -wa -wb -a sample_3_sorted.bed3+ -b sample_5_sorted.bed3+ >intersect_3_5.tsv
bedtools intersect -wa -wb -a sample_4_sorted.bed3+ -b sample_5_sorted.bed3+ >intersect_4_5.tsv


# 11. snpEff annotation
mkdir snpEff_annotation && cd snpEff_annotation
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/{print "chr"$1,$2,"rs_"NR,$3,$4,".","PASS",$5}' ../test.assignment.txt >snp.vcf
sed -i "1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" snp.vcf
sed -i 's/-/./;s/chrMT/chrM/' snp.vcf
#wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip
#unzip snpEff_v4_3_GRCh37.75.zip
#mv data ../tools/snpEff/
cd ../tools/snpEff
mkdir -p data/hg19.my && cd data/hg19.my
ln -s ../../../../gene/refGene.txt genes.refseq
ln -s /mnt/share/share/data/fna/hg19/all.fa sequences.fa
echo "# my genome, hg19.my" >>../../snpEff.config
echo "hg19.my.genome : Human" >>../../snpEff.config
cd -
/mnt/share/share/tool/jdk/bin/java -Xmx4g -jar ../tools/snpEff/snpEff.jar build -c ../tools/snpEff/snpEff.config -refseq -v hg19.my > snpEff_build.log 2> snpEff_build.err
/mnt/share/share/tool/jdk/bin/java -Xmx4g -jar ../tools/snpEff/snpEff.jar -c ../tools/snpEff/snpEff.config -v -canon hg19.my snp.vcf >snpEff.vcf 2>snpEff.err


# 12. BMR
mkdir bmr && cd bmr
# 1) exon, intron region
#awk 'BEGIN{FS=OFS="\t"} $2~/NM/' ../gc/refGene_longest.gpe >refGene_NM.gpe
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i ../gc/refGene_longest.gpe >intron.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i refGene_NM.gpe >intron.bed
sort -k1,1 -k2,2n intron.bed >intron_sorted.bed
rm intron.bed
bedtools merge -i intron_sorted.bed >intron_merge.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e ../gc/refGene_longest.gpe >exon.bed
#gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e refGene_NM.gpe >exon.bed
sort -k1,1 -k2,2n exon.bed >exon_sorted.bed
rm exon.bed
bedtools merge -i exon_sorted.bed >exon_merge.bed
bedtools subtract -a intron_merge.bed -b exon_merge.bed >intron_substract_exon.bed
# 2) repeat region
zcat ../repeat_masker/rmsk.txt.gz | awk 'BEGIN{FS=OFS="\t"} $6!~/_/ {print $6,$7,$8,"name","0",$10;}' >repeat.bed6
sort -k1,1 -k2,2n repeat.bed6 >repeat_sorted.bed6
bedtools merge -i repeat_sorted.bed6 >repeat_merge.bed3
rm repeat.bed6
# 3) (repeat, non-repeat) && (exon, intron)
bedtools intersect -wb -a repeat_merge.bed3 -b exon_merge.bed | cut -f1-3 >repeat_exon.bed
bedtools intersect -wb -a repeat_merge.bed3 -b intron_substract_exon.bed | cut -f1-3 >repeat_intron.bed
grep -v "_" /mnt/share/share/data/chr.size/hg19.size >hg19.size
bedtools complement -i repeat_merge.bed3 -g hg19.size >nonrepeat.bed
bedtools intersect -wb -a nonrepeat.bed -b exon_merge.bed | cut -f1-3 >nonrepeat_exon.bed
bedtools intersect -wb -a nonrepeat.bed -b intron_substract_exon.bed | cut -f1-3 >nonrepeat_intron.bed
# 4) GC and TC position in genome
../bin/motif_finder.pl /mnt/share/share/data/fna/hg19/all.fa > CG.bed3
sed -i "/_/d" CG.bed3
../bin/motif_finder.pl -m TC /mnt/share/share/data/fna/hg19/all.fa > TC.bed3
sed -i "/_/d" TC.bed3
# 5) GC && (repeat, non-repeat) && (exon, intron)
ls ../gc/region_gc_sorted_0* | while read file
  do
    sample=$( basename $file .bed3+ | sed 's/region_gc_sorted_//' )
    for type in repeat nonrepeat
      do
	bedtools intersect -wb -a $file -b "$type"_exon.bed | cut -f1-4 >gc_"$type"_"$sample"_exon.bed
	bedtools intersect -wb -a gc_"$type"_"$sample"_exon.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_exon_CG.bed
	bedtools intersect -wb -a gc_"$type"_"$sample"_exon.bed -b TC.bed3 | cut -f1-4 >gc_"$type"_"$sample"_exon_TC.bed
	bedtools subtract -a gc_"$type"_"$sample"_exon.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_exon_nonCG.bed
	bedtools subtract -a gc_"$type"_"$sample"_exon_nonCG.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_exon_nonCG_nonTC.bed
	bedtools intersect -wb -a $file -b "$type"_intron.bed | cut -f1-4 >gc_"$type"_"$sample"_intron.bed
	bedtools intersect -wb -a gc_"$type"_"$sample"_intron.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_intron_CG.bed
	bedtools intersect -wb -a gc_"$type"_"$sample"_intron.bed -b TC.bed3 | cut -f1-4 >gc_"$type"_"$sample"_intron_TC.bed
	bedtools subtract -a gc_"$type"_"$sample"_intron.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_intron_nonCG.bed
	bedtools subtract -a gc_"$type"_"$sample"_intron_nonCG.bed -b CG.bed3 | cut -f1-4 >gc_"$type"_"$sample"_intron_nonCG_nonTC.bed
      done
  done
# 6) snp without indel
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/ && $3!="-" && $4!="-"{len1=length($3); len2=length($4); if(len1==1 && len2==1) print "chr"$1,$2-1,$2;}' ../test.assignment.txt >snp.bed3
sort -k1,1 -k2,2n snp.bed3 >snp_sorted.bed3
rm snp.bed3 
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/ && $3!="-" && $4!="-"{len1=length($3); len2=length($4); if(len1==1 && len2==1) print "chr"$1,$2-1,$2,$3;}' ../test.assignment.txt >snp_ref.bed3+
sort -k1,1 -k2,2n snp_ref.bed3+ >snp_ref_sorted.bed3+
rm snp_ref.bed3+
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/ && $3=="-" || $4=="-"{print "chr"$1,$2-1,$2,$3;}' ../test.assignment.txt >indel.bed3+
awk 'BEGIN{FS=OFS="\t"} $1!~/GL/ && $3!="-" && $4!="-"{len1=length($3); len2=length($4); if(len1>1 || len2>1) print "chr"$1,$2-1,$2,$3;}' ../test.assignment.txt >>indel.bed3+
# 7) overlap with snp 
ls *_CG.bed *TC.bed | while read file
  do
    sample=$( echo $file | sed "s/gc_//" )
    bedtools intersect -wb -a $file -b snp_ref_sorted.bed3+ | cut -f1-4,8 >intersect_"$sample"
  done
# 8) mutation count for each base
ls intersect_* | while read file
  do
    echo $file
    awk 'BEGIN{FS=OFS="\t"} {if($5=="A"){sumA+=1;} else if($5=="C"){sumC+=1;} else if($5=="G"){sumG+=1;} else if($5=="T"){sumT+=1;}} END{print "A="sumA,"C="sumC,"G="sumG,"T="sumT;}' $file
  done
ls gc*_CG.bed *TC.bed | while read file
  do
    echo $file
    bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed $file | awk 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$7;sumC+=$8;sumG+=$9;sumT+=$10;} END{print "A="sumA,"C="sumC,"G="sumG,"T="sumT;}'
  done
# 9) synonymous mutation
grep synonymous ../snpEff_annotation/snpEff.vcf | awk 'BEGIN{FS=OFS="\t"} {len1=length($4); len2=length($5); if(len1==1 && len2==1) {print $1,$2-1,$2,$4,$5;}}' >snp_synonymous.bed
ls gc*_CG.bed gc*TC.bed | grep exon | while read file
  do
    sample=$( echo $file | sed "s/gc_//" )
    bedtools intersect -wb -a $file -b snp_synonymous.bed | cut -f1-4,8 >intersect_sysnonymous_"$sample"
  done
ls intersect_sysnonymous_* | while read file
  do
    echo $file
    awk 'BEGIN{FS=OFS="\t"} {if($5=="A"){sumA+=1;} else if($5=="C"){sumC+=1;} else if($5=="G"){sumG+=1;} else if($5=="T"){sumT+=1;}} END{print "A="sumA,"C="sumC,"G="sumG,"T="sumT;}' $file
  done


# 13. driver gene
mkdir driver_gene && cd driver_gene
# 1) CG && (repeat, nonrepeat) && (exon, intron)
bedtools intersect -a ../bmr/CG.bed3 -b ../bmr/repeat_merge.bed3 >CG_repeat.bed
bedtools intersect -a CG_repeat.bed -b ../bmr/exon_merge.bed >CG_repeat_exon.bed
bedtools intersect -a CG_repeat.bed -b ../bmr/intron_substract_exon.bed >CG_repeat_intron.bed
bedtools subtract -a ../bmr/CG.bed3 -b ../bmr/repeat_merge.bed3 > CG_nonrepeat.bed
bedtools intersect -a CG_nonrepeat.bed -b ../bmr/exon_merge.bed >CG_nonrepeat_exon.bed
bedtools intersect -a CG_nonrepeat.bed -b ../bmr/intron_substract_exon.bed >CG_nonrepeat_intron.bed
# 2) TC && (repeat, nonrepeat) && (exon, intron)
bedtools intersect -a ../bmr/TC.bed3 -b ../bmr/repeat_merge.bed3 >TC_repeat.bed
bedtools intersect -a TC_repeat.bed -b ../bmr/exon_merge.bed >TC_repeat_exon.bed
bedtools intersect -a TC_repeat.bed -b ../bmr/intron_substract_exon.bed >TC_repeat_intron.bed
bedtools subtract -a ../bmr/TC.bed3 -b ../bmr/repeat_merge.bed3 > TC_nonrepeat.bed
bedtools intersect -a TC_nonrepeat.bed -b ../bmr/exon_merge.bed >TC_nonrepeat_exon.bed
bedtools intersect -a TC_nonrepeat.bed -b ../bmr/intron_substract_exon.bed >TC_nonrepeat_intron.bed
# 3) non_CG_nonTC && (repeat, nonrepeat) && (exon, intron)
bedtools subtract -a ../bmr/exon_merge.bed -b ../bmr/CG.bed3 >exon_nonCG.bed
bedtools subtract -a exon_nonCG.bed -b ../bmr/TC.bed3 >exon_nonCG_nonTC.bed
bedtools intersect -a exon_nonCG_nonTC.bed -b ../bmr/repeat_merge.bed3 >exon_nonCG_nonTC_repeat.bed
bedtools intersect -a exon_nonCG_nonTC.bed -b ../bmr/nonrepeat.bed >exon_nonCG_nonTC_nonrepeat.bed
bedtools subtract -a ../bmr/intron_substract_exon.bed -b ../bmr/CG.bed3 >intron_nonCG.bed
bedtools subtract -a intron_nonCG.bed -b ../bmr/TC.bed3 >intron_nonCG_nonTC.bed
bedtools intersect -a intron_nonCG_nonTC.bed -b ../bmr/repeat_merge.bed3 > intron_nonCG_nonTC_repeat.bed
bedtools intersect -a intron_nonCG_nonTC.bed -b ../bmr/nonrepeat.bed > intron_nonCG_nonTC_nonrepeat.bed
# 4) split
mkdir CG
for i in {1..22} X Y M
  do
    for file in CG_repeat_exon.bed CG_repeat_intron.bed CG_nonrepeat_exon.bed CG_nonrepeat_intron.bed
      do
        grep -w "^chr"$i $file >CG/"chr""$i"_"$file"
      done
  done
mkdir TC
for i in {1..22} X Y M
  do
    for file in TC_repeat_exon.bed TC_repeat_intron.bed TC_nonrepeat_exon.bed TC_nonrepeat_intron.bed
      do
        grep -w "^chr"$i $file >TC/"chr""$i"_"$file"
      done
  done
mkdir nonCG_nonTC
for i in {1..22} X Y M
  do
    for file in exon_nonCG_nonTC_repeat.bed exon_nonCG_nonTC_nonrepeat.bed intron_nonCG_nonTC_repeat.bed intron_nonCG_nonTC_nonrepeat.bed
      do
        grep -w "^chr"$i $file >nonCG_nonTC/"chr""$i"_"$file"
      done
  done
# 5) probability
rm sumPro.tsv
cat ../gc/refGene_longest.gpe | while read line
  do
    tranName=$( echo "$line" | cut -f2 )
    geneName=$( echo "$line" | cut -f13 )
    chr=$( echo "$line" | cut -f3 )
    gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e <( echo "$line" ) >case_exon.bed
    gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i <( echo "$line" ) >case_intron.bed
    for type_1 in exon intron
      do
        if [ ! -s case_"$type_1".bed ]; then
	  for type_2 in repeat nonrepeat
	    do
	      echo "" >case_CG_"$type_2"_"$type_1".bed
	      echo "" >intersect_CG_"$type_2"_"$type_1".bed
	      echo "" >case_TC_"$type_2"_"$type_1".bed
	      echo "" >intersect_TC_"$type_2"_"$type_1".bed
	      echo "" >case_nonCG_nonTC_"$type_2"_"$type_1".bed
	      echo "" >intersect_nonCG_nonTC_"$type_2"_"$type_1".bed
	    done
	else
	  for type_2 in repeat nonrepeat
	    do
	      bedtools intersect -a case_"$type_1".bed -b CG/"$chr"_CG_"$type_2"_"$type_1".bed >case_CG_"$type_2"_"$type_1".bed
	      if [ ! -s case_CG_"$type_2"_"$type_1".bed ]; then
	        echo "" >intersect_CG_"$type_2"_"$type_1".bed
		echo "" >case_CG_"$type_2"_"$type_1".bed
	      else
	        bedtools intersect -wb -a case_CG_"$type_2"_"$type_1".bed -b ../bmr/snp_ref_sorted.bed3+ >intersect_CG_"$type_2"_"$type_1".bed
	      fi 
	      bedtools intersect -a case_"$type_1".bed -b TC/"$chr"_TC_"$type_2"_"$type_1".bed >case_TC_"$type_2"_"$type_1".bed
	      if [ ! -s case_TC_"$type_2"_"$type_1".bed ]; then
	        echo "" >intersect_TC_"$type_2"_"$type_1".bed
		echo "" >case_TC_"$type_2"_"$type_1".bed
	      else
	        bedtools intersect -wb -a case_TC_"$type_2"_"$type_1".bed -b ../bmr/snp_ref_sorted.bed3+ >intersect_TC_"$type_2"_"$type_1".bed
	      fi
	      bedtools intersect -a case_"$type_1".bed -b nonCG_nonTC/"$chr"_"$type_1"_nonCG_nonTC_"$type_2".bed >case_nonCG_nonTC_"$type_2"_"$type_1".bed
	      if [ ! -s case_nonCG_nonTC_"$type_2"_"$type_1".bed ]; then
	        echo "" >intersect_nonCG_nonTC_"$type_2"_"$type_1".bed
		echo "" >case_nonCG_nonTC_"$type_2"_"$type_1".bed
	      else
	        bedtools intersect -wb -a case_nonCG_nonTC_"$type_2"_"$type_1".bed -b ../bmr/snp_ref_sorted.bed3+ >intersect_nonCG_nonTC_"$type_2"_"$type_1".bed
	      fi
	    done
	fi
      done
    exon_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_nonCG_nonTC_repeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    exon_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_nonCG_nonTC_nonrepeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    intron_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_nonCG_nonTC_repeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    intron_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_nonCG_nonTC_nonrepeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    cg_exon_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_CG_repeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    cg_exon_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_CG_nonrepeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    cg_intron_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_CG_repeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    cg_intron_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_CG_nonrepeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    tc_exon_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_TC_repeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    tc_exon_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_TC_nonrepeat_exon.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    tc_intron_repeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_TC_repeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    tc_intron_nonrepeat_bases=$( bedtools nuc -fi /mnt/share/share/data/fna/hg19/all.fa -bed case_TC_nonrepeat_intron.bed | awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} NR>1{sumA+=$11; sumC+=$12; sumG+=$13; sumT+=$14;} END{print sumA,sumC,sumG,sumT}' ) 
    mut_exon_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_nonCG_nonTC_repeat_exon.bed )
    mut_exon_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_nonCG_nonTC_nonrepeat_exon.bed )
    mut_intron_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_nonCG_nonTC_repeat_intron.bed )
    mut_intron_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_nonCG_nonTC_nonrepeat_intron.bed )
    mut_cg_exon_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_CG_repeat_exon.bed )
    mut_cg_exon_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_CG_nonrepeat_exon.bed )
    mut_cg_intron_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_CG_repeat_intron.bed )
    mut_cg_intron_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_CG_nonrepeat_intron.bed )
    mut_tc_exon_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_TC_repeat_exon.bed )
    mut_tc_exon_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_TC_nonrepeat_exon.bed )
    mut_tc_intron_repeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_TC_repeat_intron.bed )
    mut_tc_intron_nonrepeat_bases=$( awk -v sumA=0 -v sumC=0 -v sumG=0 -v sumT=0 'BEGIN{FS=OFS="\t"} { if($12=="A"){sumA+=1;} else if($12=="C"){sumC+=1;} else if($12=="G"){sumG+=1;} else if($12=="T"){sumT+=1;}} END{print sumA,sumC,sumG,sumT}' intersect_TC_nonrepeat_intron.bed )
     ./binomTest.R --exon_repeat_bases "$exon_repeat_bases" --exon_nonrepeat_bases "$exon_nonrepeat_bases" --intron_repeat_bases "$intron_repeat_bases" --intron_nonrepeat_bases "$intron_nonrepeat_bases" --cg_exon_repeat_bases "$cg_exon_repeat_bases" --cg_exon_nonrepeat_bases "$cg_exon_nonrepeat_bases" --cg_intron_repeat_bases "$cg_intron_repeat_bases" --cg_intron_nonrepeat_bases "$cg_intron_nonrepeat_bases" --tc_exon_repeat_bases "$tc_exon_repeat_bases" --tc_exon_nonrepeat_bases "$tc_exon_nonrepeat_bases" --tc_intron_repeat_bases "$tc_intron_repeat_bases" --tc_intron_nonrepeat_bases "$tc_intron_nonrepeat_bases" --mut_exon_repeat_bases "$mut_exon_repeat_bases" --mut_exon_nonrepeat_bases "$mut_exon_nonrepeat_bases" --mut_intron_repeat_bases "$mut_intron_repeat_bases" --mut_intron_nonrepeat_bases "$mut_intron_nonrepeat_bases" --mut_cg_exon_repeat_bases "$mut_cg_exon_repeat_bases" --mut_cg_exon_nonrepeat_bases "$mut_cg_exon_nonrepeat_bases" --mut_cg_intron_repeat_bases "$mut_cg_intron_repeat_bases" --mut_cg_intron_nonrepeat_bases "$mut_cg_intron_nonrepeat_bases" --mut_tc_exon_repeat_bases "$mut_tc_exon_repeat_bases" --mut_tc_exon_nonrepeat_bases "$mut_tc_exon_nonrepeat_bases" --mut_tc_intron_repeat_bases "$mut_tc_intron_repeat_bases" --mut_tc_intron_nonrepeat_bases "$mut_tc_intron_nonrepeat_bases" --tranName "$tranName" --geneName "$geneName" >>sumPro.tsv
   done
<<mark
sort -k3,3n sumPro.tsv >sumPro_sorted.tsv
rm sumPro.tsv
awk -v myIndex=0 -v myCom="NA" 'BEGIN{FS=OFS="\t"} NR==1{myIndex=1; myCom=$3; print $0,myIndex;} NR>1{if($3!=myCom){myIndex+=1; myCom=$3;}  print $0,myIndex;}' sumPro_sorted.tsv >sumPro_rank.tsv
[ -e score.tsv ] && rm score.tsv
cat sumPro_rank.tsv | while read line
  do
    myPro=$( echo "$line" | cut -f3 )
    myRank=$( echo "$line" | cut -f4 )
    myScore=$( echo "-(l(28108*$myPro/$myRank)/l(10))" | bc -l )
    echo -e "$line\t$myScore" >>score.tsv
  done
mark
# 6) fisher's combind p-value
./pchisq.R sumPro.tsv >pchisq.tsv
# 7) fdr
./fdr.R pchisq.tsv >fdr.tsv
# 8) driver gene
awk '$4<0.05' fdr.tsv | sort -k4,4n >driver_gene.tsv
#selectionOutput.pl --colmun12 2 driver_gene.tsv refGene.gpe - 2>/dev/null | awk 'BEGIN{FS=OFS="\t"} {print $3,$5,$6,$2;}' >driver_refGene.bed
[ -e driver_vcf.tsv ] && rm driver_vcf.tsv
cat driver_gene.tsv | while read line
  do
    tranName=$( echo "$line" | cut -f1 )
    grep -w $tranName ../snpEff_annotation/snpEff.vcf | while read vcf
      do
	echo -e "$line\t$vcf" >>driver_vcf.tsv
      done
  done

# 9) snpEff
perl -e '
	open IN, $ARGV[0];
	while(<IN>){
		chomp;
		my @input=split "\t",$_;
		my $snpEff=$input[11];
		my @SNPEffs=split ",",$snpEff;
		my $effCount=scalar(@SNPEffs);
		if($effCount==1){
			print $_."\n";
		}else{
			my $tranName=$input[0];
			my $output=join "\t",@input[0..10];
			for(my $i=0; $i<$effCount; $i++){
				my @detailFun=split "\\|", $SNPEffs[$i];
				my $comName=$detailFun[6];
				$comName=~s/\..*//;
				next if $comName ne $tranName;
				$output.="\t".$SNPEffs[$i];
			}
			print $output."\n";
		}
	}
' driver_vcf.tsv >driver_vcf_process.tsv
cut -f12 driver_vcf_process.tsv | cut -d "|" -f2 | sort | uniq -c


# 10) rpkm
wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2refseq.gz
gunzip gene2refseq.gz
grep "^9606" gene2refseq >human_gene2refseq.tsv
awk 'BEGIN{FS=OFS="\t"} {sub("\\..*","",$4); print $2,$4,$16;}' human_gene2refseq.tsv >human_gene2refseq_process.tsv
sort human_gene2refseq_process.tsv | uniq >human_gene2refseq_process_1.tsv
mv human_gene2refseq_process_1.tsv human_gene2refseq_process.tsv
rm human_gene2refseq.tsv
selectionOutput.pl --colmun12 10 driver_gene.tsv RPKM.tsv >driver_rpkm.tsv 2>/dev/null
grep -i colon driver_rpkm.tsv >driver_rpkm_colon.tsv
selectionOutput.pl --colmun11 10 driver_rpkm_colon.tsv driver_gene.tsv >/dev/null 2>driver_gene_noRPKM.tsv
selectionOutput.pl --colmun11 2 --colmun12 2 driver_gene_noPRKM.tsv RPKM.tsv >driver_rpkm_symbol.tsv 2>/dev/null
selectionOutput.pl --colmun11 2 --colmun12 2 driver_gene_noRPKM.tsv RPKM.tsv >driver_rpkm_symbol.tsv 2>/dev/null
grep -i colon driver_rpkm_symbol.tsv >driver_rpkm_symbol_colon.tsv
rm driver_rpkm_symbol.tsv
perl -e '
	open IN_1, $ARGV[0];
	open IN_2, $ARGV[1];
	my %hash;
	while(<IN_1>){
		chomp;
		my ($tranName, $geneName)=(split "\t", $_)[0,1];
		$hash{$geneName}=$tranName;
	}
	while(<IN_2>){
		chomp;
		my @input=split "\t", $_;
		my $geneName=$input[1];
		$input[9]=$hash{$geneName};
		print join "\t",@input;
		print "\n";
	}
' driver_gene_noRPKM.tsv driver_rpkm_symbol_colon.tsv >driver_rpkm_symbol_colon_replace.tsv
cat driver_rpkm_symbol_colon_replace.tsv >>driver_rpkm_colon.tsv
selectionOutput.pl --colmun11 10 driver_rpkm_colon.tsv driver_gene.tsv >/dev/null 2>driver_gene_noRPKM.tsv
[ -e driver_gene_new_trans.tsv ] && rm driver_gene_new_trans.tsv
cat driver_gene_noRPKM.tsv | while read line_1
  do
    tranName=$( echo "$line_1" | cut -f1 )
    geneName=$( echo "$line_1" | cut -f2 )
    grep -w $geneName human_gene2refseq_process.tsv | while read line_2
      do
        tranName_2=$( echo "$line_2" | cut -f2 )
	grepCount=$( grep -w "$tranName_2" RPKM.tsv | head -n 1 | wc -l )
        if [ "$grepCount" -eq "1" ]; then
	  echo -e "$line_1\t$line_2" >>driver_gene_new_trans.tsv
	fi
      done
  done
perl -e '
	open IN_1, $ARGV[0];
	open IN_2, $ARGV[1];
	my %hash;
	while(<IN_1>){
		chomp;
		my ($tranName_1, $tranName_2)=(split "\t", $_)[0,5];
		$hash{$tranName_2}=$tranName_1;
	}
	while(<IN_2>){
		chomp;
		my @input=split "\t", $_;
		my $tranName=$input[9];
		if(exists $hash{$tranName}){
			$input[9]=$hash{$tranName};
			print join "\t",@input;
			print "\n";
		}
	}
' driver_gene_new_trans.tsv RPKM.tsv >driver_rpkm_new_trans.tsv
grep -i colon driver_rpkm_new_trans.tsv >>driver_rpkm_colon.tsv
awk '$4>1' driver_rpkm_colon.tsv | cut -f10 | sort | uniq | wc -l

# 2) expression gene with high mutaiton sites
awk '$4>1' driver_rpkm_colon.tsv | cut -f10 >expression_colon_gene.tsv
selectionOutput.pl expression_colon_gene.tsv driver_vcf_process.tsv >expression_colon_trans.tsv 2>/dev/null
grep -P "stop_gained|stop_lost" expression_colon_trans.tsv | cut -f2 | sort | uniq | wc -l


# 11) cancer gene
cp ../gene/cancer_gene_census.csv.bak cancer_gene_census.csv
cut -d "," -f1 cancer_gene_census.csv >cancer_gene_name.tsv
sed -i '1d' cancer_gene_name.tsv 
selectionOutput.pl --colmun12 13 cancer_gene_name.tsv refGene.gpe >cancer_gene_refseq.tsv 2>/dev/nul
mkdir cancer_gene && cd cancer_gene
ln -s ../cancer_gene_refseq.tsv part.gpe
rm sumPro.tsv
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e part.gpe >case_exon_cancer.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i part.gpe >case_intron_cancer.bed
./binomTest_cancer_gene.sh >binomTest_cancer_gene.log 2>binomTest_cancer_gene.err

selectionOutput.pl --colmun11 2 part.gpe ../driver_gene.tsv >cancer_driver.tsv 2>/dev/null


# 12) methylation
mkdir methy && cd methy
bedtools intersect -a ../case_exon_cancer.bed -b ../../../HAIB_Methyl/methy_merge_sorted.bed3 >case_exon_methy.bed
bedtools intersect -a ../case_intron_cancer.bed -b ../../../HAIB_Methyl/methy_merge_sorted.bed3 >case_intron_methy.bed
./binomTest_cancer_methy.sh >binomTest_cancer_methy.log 2>binomTest_cancer_methy.err
mkdir all_genes
mv case* intersect* all_gene
selectionOutput.pl --colmun12 2 ../cancer_driver.tsv ../part.gpe >part.gpe 2>/dev/null
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e part.gpe >case_exon_cancer.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i part.gpe >case_intron_cancer.bed
bedtools intersect -a case_exon_cancer.bed -b ../../../HAIB_Methyl/methy_merge_sorted.bed3 >case_exon_methy.bed
bedtools intersect -a case_intron_cancer.bed -b ../../../HAIB_Methyl/methy_merge_sorted.bed3 >case_intron_methy.bed
./binomTest_cancer_methy.sh >binomTest_cancer_methy.log 2>binomTest_cancer_methy.err

# 13) enhancer
mkdir enhancer && cd enhancer
bedtools intersect -a ../case_exon_cancer.bed -b ../../../enhancer/enhancer_merge_sorted.bed3 >case_exon_methy.bed
bedtools intersect -a ../case_intron_cancer.bed -b ../../../enhancer/enhancer_merge_sorted.bed3 >case_intron_methy.bed
./binomTest_cancer_methy.sh >binomTest_cancer_methy.log 2>binomTest_cancer_methy.err
mkdir all_genes
mv case* intersect* all_gene
ln -s ../methy/part.gpe .
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -e part.gpe >case_exon_cancer.bed
gpeFeature.pl -b -g /mnt/share/share/data/chr.size/hg19.size -i part.gpe >case_intron_cancer.bed
bedtools intersect -a case_exon_cancer.bed -b ../../../enhancer/enhancer_merge_sorted.bed3 >case_exon_methy.bed
bedtools intersect -a case_intron_cancer.bed -b ../../../enhancer/enhancer_merge_sorted.bed3 >case_intron_methy.bed
./binomTest_cancer_methy.sh >binomTest_cancer_methy.log 2>binomTest_cancer_methy.err
