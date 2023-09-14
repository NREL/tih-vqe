#!/bin/bash

# finish runs that didn't have time to converge
ucc=( s sd )

# array=( 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 )
array=( 10 )

for u in "${ucc[@]}"; do
  for i in "${array[@]}"; do
    mv parallel.log-$u-$i parallel.log-$u-$i\_old
    sbatch submit.sh $i $u
  done
done
