#!/bin/bash

usage(){
 cat << EOF
Description:
    This script was used to calculated the time between two files.
Usage:
    timeCalculated.sh file_1 file_2
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

file_1=$( readlink -f $1 )
file_2=$( readlink -f $2 )
time_1=$( stat -c %Y $file_1 )
time_2=$( stat -c %Y $file_2 )
timeCal=$(( $time_1 - $time_2 ))
timeCal=${timeCal#-}

timeCalFormat=$(show_time $timeCal)

echo -e "${timeCal}\t${timeCalFormat}"
