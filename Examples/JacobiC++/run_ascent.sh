#!/bin/bash -l
#SBATCH --job-name="ascent"
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --constraint=gpu
#SBATCH --account=csstaff
##SBATCH --reservation=insitu
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --partition=debug
##SBATCH --mail-type=ALL
##SBATCH --mail-user=jfavre@cscs.ch
##SBATCH --output=/users/jfavre/out.log
##SBATCH --error=/users/jfavre/err.log

module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel
module use /scratch/snx3000tds/jfavre/daint/modules/all

module load Ascent/0.9.1-CrayGNU-21.09

mkdir -p $SCRATCH/Ascent/test
pushd    $SCRATCH/Ascent/test

cp /users/jfavre/Projects/InSitu/InSitu-Vis-Tutorial2022/Examples/JacobiC++/{ascent_actions.yaml,plot_actions.yaml} $SCRATCH/Ascent/test

srun -n $SLURM_NNODES  /users/jfavre/Projects/InSitu/InSitu-Vis-Tutorial2022/Examples/JacobiC++/buildAscent/bin/pjacobi  --res=256 --mesh=uniform
