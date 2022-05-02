#!/bin/bash -l
#SBATCH --job-name="sedov-catalyst"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
#SBATCH --account=csstaff
#SBATCH --partition=debug
#SBATCH --exclusive
#SBATCH --time=00:10:00
##SBATCH --mail-type=ALL
##SBATCH --mail-user=jfavre@cscs.ch
##SBATCH --output=/users/jfavre/out.log
##SBATCH --error=/users/jfavre/err.log
#SBATCH --hint=nomultithread


export OMP_NUM_THREADS=12
export OMP_PLACES=sockets
export OMP_PROC_BIND=close

module load daint-gpu
module load ParaView

mkdir -p $SCRATCH/SPH-EXA/Catalyst
pushd $SCRATCH/SPH-EXA/Catalyst

srun -n $SLURM_NNODES /users/jfavre/Projects/SPH-EXA/buildCatalystDaint/main/src/sphexa/sphexa  --init sedov \
                      --quiet -s 100 -n 100 \
                      --catalyst /users/jfavre/Projects/SPH-EXA/runs/catalyst_helloworld.py
