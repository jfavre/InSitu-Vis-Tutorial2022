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

module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu
module load Ascent/0.9.1-CrayGNU-21.09

mkdir -p $SCRATCH/Ascent/test

cwd=$PWD
cp $PWD/{ascent_actions.yaml,plot_actions.yaml} $SCRATCH/Ascent/test

pushd    $SCRATCH/Ascent/test
srun -n $SLURM_NNODES  $cwd/buildAscent/bin/pjacobi  --res=128 --mesh=uniform

popd
cp $SCRATCH/Ascent/test/Jacobi.* .
