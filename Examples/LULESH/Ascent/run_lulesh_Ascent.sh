#!/bin/bash -l
#SBATCH --job-name="lulesh+ascent"
#SBATCH --account="csstaff"
##SBATCH --reservation="in-situ"
#SBATCH --time=00:10:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

# should run one task per node becasue Asent needs a dedicated CUDA node per task

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir $SCRATCH/Lulesh

cp ../buildAscent/bin/lulesh2.0 $SCRATCH/Lulesh

# enable a custom scene description
# cp ascent_actions.yaml trigger_ascent_actions.yaml  $SCRATCH/Lulesh

pushd $SCRATCH/Lulesh
srun ./lulesh2.0 -s 30 -p -i 1000
