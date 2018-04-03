#!/bin/bash

if [ -z sumPro.tsv ]; then
  rm sumPro.tsv
fi

#cat ../gc/refGene_longest.gpe | while read line
cat test.gpe | while read line
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
