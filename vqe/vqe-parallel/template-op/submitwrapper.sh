#!/bin/bash

# change which op to run in submit.sh
ucc=( s sd )

# bond length index array
array=( 30 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 )

for u in "${ucc[@]}"; do
  for i in "${array[@]}"; do
    sbatch submit.sh $i $u
  done
done
