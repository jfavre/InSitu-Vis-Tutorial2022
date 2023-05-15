#!/bin/bash -l
#SBATCH --job-name="lulesh+ascent+paraview"
##SBATCH --account="csstaff"
##SBATCH --reservation="in-situ"
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

# should run one task per node because Ascent needs a dedicated CUDA node per task
module load cudatoolkit/21.3_11.2 Ascent/0.9.1-CrayGNU-21.09

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p $SCRATCH/Lulesh/Ascent

cp ../../buildAscent/bin/lulesh2.0 $SCRATCH/Lulesh/Ascent

# enable a custom scene description
cp ascent_actions.yaml paraview-pressure-threshold.py  $SCRATCH/Lulesh/Ascent
#cp trigger_isocontours_actions.yaml  $SCRATCH/Lulesh

pushd $SCRATCH/Lulesh/Ascent
srun ./lulesh2.0 -s 30 -p -i 1000
