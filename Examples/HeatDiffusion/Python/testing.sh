#!/bin/bash
#SBATCH --job-name="Ascent-Python-Jacobi"
#SBATCH --nodes=1
#SBATCH --ntasks=2

#SBATCH --account=usup
#SBATCH --time=00:05:00
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --exclusive

module load daint-gpu
module load Ascent
module load matplotlib

mkdir $SCRATCH/TAscent
cp jacobi_insitu_serial.py jacobi_insitu_Ascent.py jacobi_insitu_parallel.py jacobi_insitu_parallel_Ascent.py $SCRATCH/TAscent

pushd $SCRATCH/TAscent
srun -n 1 python3 jacobi_insitu_serial.py
srun -n 1 python3 jacobi_insitu_Ascent.py
srun -n $SLURM_NTASKS -N $SLURM_NNODES python3 jacobi_insitu_parallel.py
srun -n $SLURM_NTASKS -N $SLURM_NNODES python3 jacobi_insitu_parallel_Ascent.py
