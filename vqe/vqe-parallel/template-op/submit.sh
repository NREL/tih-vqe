#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=36 --partition=debug --time=1:00:00 --qos=normal
#SBATCH --account=scatter
#SBATCH --exclusive
#SBATCH --export=ALL

module load conda
source ~/conda_init
conda activate dompi
module load intel-mpi/2020.1.217

i=$1
ucc=$2
if [ -z $i ]; then
    echo "i defaulting to 1"
    i=1
else echo "i is set to '$i'"
fi

rm parallel.log-$i

op=2 # change which op to run in vqeops.py
./vqe-H-copy.py $op $SLURM_NTASKS $i

srun ./vqe-wrapper.py $op $i $ucc > parallel.log-$ucc-$i
