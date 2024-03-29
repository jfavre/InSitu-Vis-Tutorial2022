#!/bin/bash -l
#SBATCH --job-name="lulesh+catalyst"
##SBATCH --account="csstaff"
##SBATCH --reservation="insitu"
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu
module load ParaView/5.11.1-CrayGNU-21.09-EGL
export CATALYST_IMPLEMENTATION_PATHS=$EBROOTPARAVIEW/lib64/catalyst

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export VTK_SILENCE_GET_VOID_POINTER_WARNINGS=1

mkdir -p $SCRATCH/Lulesh/Catalyst

cp catalyst.py ../buildCatalyst/bin/lulesh2.0 $SCRATCH/Lulesh/Catalyst

pushd $SCRATCH/Lulesh/Catalyst

srun ./lulesh2.0 -x catalyst.py -s 30 -p -i 1000

