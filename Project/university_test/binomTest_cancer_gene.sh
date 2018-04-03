#!/bin/bash

[ -e sumPro.tsv ] && rm sumPro.tsv
rm case_CG_*
rm case_TC_*
rm case_nonCG_*

tranName="NA"
geneName="NA"

for i in {1..22} X Y
  do
    grep -w ^chr"$i" case_exon_cancer.bed >case_exon.bed
    grep -w ^chr"$i" case_intron_cancer.bed >case_intron.bed
    chr="chr"$i
    for type_1 in exon intron
      do
	for type_2 in repeat nonrepeat
	  do
	    bedtools intersect -a case_"$type_1".bed -b ../CG/"$chr"_CG_"$type_2"_"$type_1".bed >temp_case_CG_"$type_2"_"$type_1".bed
	    if [ -s temp_case_CG_"$type_2"_"$type_1".bed ]; then
	      cat temp_case_CG_"$type_2"_"$type_1".bed >>case_CG_"$type_2"_"$type_1".bed
	    fi 
	    bedtools intersect -a case_"$type_1".bed -b ../TC/"$chr"_TC_"$type_2"_"$type_1".bed >temp_case_TC_"$type_2"_"$type_1".bed
	    if [ -s case_TC_"$type_2"_"$type_1".bed ]; then
	      cat temp_case_TC_"$type_2"_"$type_1".bed >>case_TC_"$type_2"_"$type_1".bed
	    fi
	    bedtools intersect -a case_"$type_1".bed -b ../nonCG_nonTC/"$chr"_"$type_1"_nonCG_nonTC_"$type_2".bed >temp_case_nonCG_nonTC_"$type_2"_"$type_1".bed
	    if [ -s case_nonCG_nonTC_"$type_2"_"$type_1".bed ]; then
	      cat temp_case_nonCG_nonTC_"$type_2"_"$type_1".bed >>case_nonCG_nonTC_"$type_2"_"$type_1".bed
	    fi
	  done
      done
  done

for type_1 in exon intron
  do
    for type_2 in repeat nonrepeat
      do
        if [ ! -s case_CG_"$type_2"_"$type_1".bed ]; then
	  echo "" >case_CG_"$type_2"_"$type_1".bed
	  echo "" >intersect_CG_"$type_2"_"$type_1".bed
	else
	  bedtools intersect -wb -a case_CG_"$type_2"_"$type_1".bed -b ../../bmr/snp_ref_sorted.bed3+ >intersect_CG_"$type_2"_"$type_1".bed
	fi
	if [ ! -s case_TC_"$type_2"_"$type_1".bed ]; then
	  echo "">case_TC_"$type_2"_"$type_1".bed
	  echo "" >intersect_TC_"$type_2"_"$type_1".bed
	else
	  bedtools intersect -wb -a case_TC_"$type_2"_"$type_1".bed -b ../../bmr/snp_ref_sorted.bed3+ >intersect_TC_"$type_2"_"$type_1".bed
	fi
	if [ ! -s case_nonCG_nonTC_"$type_2"_"$type_1".bed ]; then
	  echo "" >case_nonCG_nonTC_"$type_2"_"$type_1".bed
	  echo "" >intersect_nonCG_nonTC_"$type_2"_"$type_1".bed
	else
	  bedtools intersect -wb -a case_nonCG_nonTC_"$type_2"_"$type_1".bed -b ../../bmr/snp_ref_sorted.bed3+ >intersect_nonCG_nonTC_"$type_2"_"$type_1".bed
        fi
      done
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
../binomTest.R --exon_repeat_bases "$exon_repeat_bases" --exon_nonrepeat_bases "$exon_nonrepeat_bases" --intron_repeat_bases "$intron_repeat_bases" --intron_nonrepeat_bases "$intron_nonrepeat_bases" --cg_exon_repeat_bases "$cg_exon_repeat_bases" --cg_exon_nonrepeat_bases "$cg_exon_nonrepeat_bases" --cg_intron_repeat_bases "$cg_intron_repeat_bases" --cg_intron_nonrepeat_bases "$cg_intron_nonrepeat_bases" --tc_exon_repeat_bases "$tc_exon_repeat_bases" --tc_exon_nonrepeat_bases "$tc_exon_nonrepeat_bases" --tc_intron_repeat_bases "$tc_intron_repeat_bases" --tc_intron_nonrepeat_bases "$tc_intron_nonrepeat_bases" --mut_exon_repeat_bases "$mut_exon_repeat_bases" --mut_exon_nonrepeat_bases "$mut_exon_nonrepeat_bases" --mut_intron_repeat_bases "$mut_intron_repeat_bases" --mut_intron_nonrepeat_bases "$mut_intron_nonrepeat_bases" --mut_cg_exon_repeat_bases "$mut_cg_exon_repeat_bases" --mut_cg_exon_nonrepeat_bases "$mut_cg_exon_nonrepeat_bases" --mut_cg_intron_repeat_bases "$mut_cg_intron_repeat_bases" --mut_cg_intron_nonrepeat_bases "$mut_cg_intron_nonrepeat_bases" --mut_tc_exon_repeat_bases "$mut_tc_exon_repeat_bases" --mut_tc_exon_nonrepeat_bases "$mut_tc_exon_nonrepeat_bases" --mut_tc_intron_repeat_bases "$mut_tc_intron_repeat_bases" --mut_tc_intron_nonrepeat_bases "$mut_tc_intron_nonrepeat_bases" --tranName "$tranName" --geneName "$geneName" >>sumPro.tsv
../pchisq.R sumPro.tsv >pchisq.tsv
