#!/bin/bash -l
#SBATCH --job-name="catalyst"
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --constraint=gpu
#SBATCH --account=csstaff
##SBATCH --reservation=insitu
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --partition=normal
##SBATCH --mail-type=ALL
##SBATCH --mail-user=jfavre@cscs.ch
##SBATCH --output=/users/jfavre/out.log
##SBATCH --error=/users/jfavre/err.log

module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu
module load ParaView/5.11.1-CrayGNU-21.09-EGL
export CATALYST_IMPLEMENTATION_PATHS=$EBROOTPARAVIEW/lib64/catalyst

mkdir -p $SCRATCH/Catalyst/test
pushd    $SCRATCH/Catalyst/test

cp $PWD/catalyst_state.py $SCRATCH/Catalyst/test

srun -n $SLURM_NNODES  $PWD/buildCatalyst/bin/pjacobi  --res=256 --mesh=uniform catalyst_state.py
