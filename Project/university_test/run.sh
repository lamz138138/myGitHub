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


# 10. miRNA
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


