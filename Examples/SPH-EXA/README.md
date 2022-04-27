SPH-EXA is a Smoothed Particle Hydrodynamics code in 3D. All details can be found at https://www.pasc-ch.org/projects/2021-2024/sph-exa2/

Source code can be found at https://github.com/unibas-dmi-hpc/SPH-EXA

the source code file ascent_adaptor.h includes a pre-defined pipeline which will draw images of all particles where rho > 1.4. o

As described in https://ascent.readthedocs.io/en/latest/Actions/Actions.html, Ascent will look in the current working directory for files called ascent_actions.json or ascent_actions.yaml. Thus, any user-defined visualization script can be defined to over-ride the pre-defined pipeline.

We provide here basic actions files.

