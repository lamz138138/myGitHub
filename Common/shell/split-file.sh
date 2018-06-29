#!/bin/bash

usage(){
 cat << EOF
Description:
    This script was used to split file to chunk.

Usage:
    split-file.sh -o split file

Options:
    -o name of directory containing output files [default: split]
    -c split file to this chunk [default: 10]

EOF
    exit 0
}

[ $1 ] || usage

outputDir=$PWD/split
chunkCount=10
while getopts "ho:c:" OPTION
 do
  case $OPTION in
        h) usage;;
        o) outputDir=$OPTARG;;
	c) chunkCount=$OPTARG;;
  esac
 done
 shift $((OPTIND - 1))

outputDir=$( readlink -f $outputDir )
inputFile=$( readlink -f $1 )
[ -e $outputDir ] && mv $outputDir ${outputDir}_bak
mkdir $outputDir && cd $outputDir

lineCount=$( wc -l $inputFile | cut -f1 -d " " )
chunk=$( echo "$lineCount/$chunkCount" | bc )
for((i=1; i<$chunkCount; i++))
  do
    start=$[ ($i-1)*$chunk+1 ]
    end=$[ $i*chunk ]
    sed -n "${start},${end}p" $inputFile >chunk-$i
  done
start=$[ ($chunkCount-1)*$chunk+1 ]
end=$lineCount
sed -n "${start},${end}p" $inputFile >chunk-$chunkCount
