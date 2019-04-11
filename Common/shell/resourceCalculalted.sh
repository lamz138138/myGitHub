#!/bin/bash

usage(){
 cat << EOF
Description:
    This script was used to calculated resource in runing pipeline (There should be no output in .err).
Usage:
    timeCalculated.sh files_to_be_listed (eg: ./*/defuse-*o*)
Options:
    -h  Help
EOF
    exit 0
}

[ $# -eq 0 ] && usage
while getopts "h" OPTION
 do
  case $OPTION in
        h) usage;;
  esac
 done
shift $((OPTIND - 1))

function show_time(){
    if [ $1 -lt 86400 ]; then 
        date -d@${1} -u '+%Hh:%Mm:%Ss';
    else 
        echo "$(($1/86400)) days $(date -d@$(($1%86400)) -u '+%Hh:%Mm:%Ss')" ;
    fi
}

[ -e temp_resource.tsv ] && rm temp_resource.tsv
for file_1 in "$@"
    do
        file_2=${file_1/\.o/\.e}
        errSize=$( ls -alh $file_2 | cut -d " " -f5 )
        if [ $errSize -ne "0" ]; then
            exit 0
        fi
        # 1. running time
        file_1=$( readlink -f $file_1 )
        file_2=$( readlink -f $file_2 )
        time_1=$( stat -c %Y $file_1 )
        time_2=$( stat -c %Y $file_2 )
        timeCal=$(( $time_1 - $time_2 ))
        timeCal=${timeCal#-}
        timeCalFormat=$(show_time $timeCal)
        # 2. cup time
        cpuT=$(grep "Resources Used" $file_1 | sed "s#.*cput=##;s#,.*##" )
        cpuS=$( echo $cpuT | awk -F'[-:]' 'NF==4{t=$4+60*($3+60*($2+24*$1));} NF==3{t=$3+60*($2+60*$1);} {print t;}' )
        cpuSFormat=$(show_time $cpuS)
        # 3. vmem
        vmem=$(grep "Resources Used" $file_1 | sed "s#.*vmem=##;s#kb,.*##" )
        # 4. mem
        mem=$(grep "Resources Used" $file_1 | sed "s#.*,mem=##;s#kb,.*##" )
        # 5. output
        echo -e "${timeCal}\t${timeCalFormat}\t${cpuS}\t${cpuSFormat}\t${vmem}\t${mem}" >>temp_resource.tsv
    done

# Get min, max, median, aveage
# a. with perl (can be extened)
# -0777 : read the whole file at once instead of line by line
# -a : autosplit into the @F array
# If you want decimals, replace %d with something like %.2f
#cut -f 1 temp_resource.tsv | perl -M'List::Util qw(sum max min)' -MPOSIX -0777 -a -ne 'printf "%-7s : %d\n"x4, "Min", min(@F), "Max", max(@F), "Average", sum(@F)/@F,  "Median", sum( (sort {$a<=>$b} @F)[ int( $#F/2 ), ceil( $#F/2 ) ] )/2;'

# b. with awk (it calculates median as mean of the two central values if value count is even)
timeCalAll=$( cut -f 1 temp_resource.tsv | sort -n | awk 'OFS="\t" {a[i++]=$0;s+=$0} END{print a[0],a[i-1],(a[int(i/2)]+a[int((i-1)/2)])/2,int(s/i+0.5);}' )
cpuSAll=$( cut -f 3 temp_resource.tsv | sort -n | awk 'OFS="\t" {a[i++]=$0;s+=$0} END{print a[0],a[i-1],(a[int(i/2)]+a[int((i-1)/2)])/2,int(s/i+0.5);}' )
vmemAll=$( cut -f 5 temp_resource.tsv | sort -n | awk 'OFS="\t" {a[i++]=$0;s+=$0} END{print a[0]"kb",a[i-1]"kb",(int((a[int(i/2)]+a[int((i-1)/2)])/2))"kb",(int(s/i+0.5))"kb";}' )
memAll=$( cut -f 6 temp_resource.tsv | sort -n | awk 'OFS="\t" {a[i++]=$0;s+=$0} END{print a[0]"kb",a[i-1]"kb",(int((a[int(i/2)]+a[int((i-1)/2)])/2))"kb",(int(s/i+0.5))"kb";}' )

timeCalMinFormat=$( echo "$timeCalAll" | cut -f1 )
timeCalMinFormat=$( show_time $timeCalMinFormat )
timeCalMaxFormat=$( echo "$timeCalAll" | cut -f2 )
timeCalMaxFormat=$( show_time $timeCalMaxFormat )
timeCalMedianFormat=$( echo "$timeCalAll" | awk 'FS="\t" {print int($3+0.5);}' )
timeCalMedianFormat=$( show_time $timeCalMedianFormat )
timeCalMeanFormat=$( echo "$timeCalAll" | awk 'FS="\t" {print int($4+0.5);}' )
timeCalMeanFormat=$( show_time $timeCalMeanFormat )

cpuSMinFormat=$( echo "$cpuSAll" | cut -f1 )
cpuSMinFormat=$( show_time $cpuSMinFormat )
cpuSMaxFormat=$( echo "$cpuSAll" | cut -f2 )
cpuSMaxFormat=$( show_time $cpuSMaxFormat )
cpuSMedianFormat=$( echo "$cpuSAll" | awk 'FS="\t" {print int($3+0.5);}' )
cpuSMedianFormat=$( show_time $cpuSMedianFormat )
cpuSMeanFormat=$( echo "$cpuSAll" | awk 'FS="\t" {print int($4+0.5);}' )
cpuSMeanFormat=$( show_time $cpuSMeanFormat )

echo -e "#\tMin\tMax\tMedian\tMean"
echo -e "real_time:\t${timeCalAll}"
echo -e "real_time_formate:\t${timeCalMinFormat}\t${timeCalMaxFormat}\t${timeCalMedianFormat}\t${timeCalMeanFormat}"
echo -e "cpuS_time:\t${cpuSAll}"
echo -e "cpuS_time_formate:\t${cpuSMinFormat}\t${cpuSMaxFormat}\t${cpuSMedianFormat}\t${cpuSMeanFormat}"
echo -e "vmem:\t$vmemAll"
echo -e "mem:\t$memAll"

echo -e "timeCal\ttimeCalFormat\tcpuS\tcupSFormat\tvmem\tmem" >temp_resource_2.tsv
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5"kb",$6"kb";}' temp_resource.tsv >temp_resource_2.tsv
mv temp_resource_2.tsv temp_resource.tsv
