mpiexec -n 2 python3 heat_diffusion_insitu_parallel_Ascent.py \
                           --res=64 --timesteps 1000 --mesh=uniform --frequency 100 --verbose
