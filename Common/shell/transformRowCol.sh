#!/bin/bash

n=1
m=$( awk 'NR==1{print NF}' $1 )
while ((n<=m))
 do
  cut -f$n $1 | tr "\n" "\t"; echo
  (( n++ ))
 done
