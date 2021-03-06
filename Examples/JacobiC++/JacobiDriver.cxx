/*
 A parallel Jacobi solver for the Laplacian equation in 2D
 Written by Jean M. Favre, Swiss National Supercomputing Center
 Last tested Tue Jul 19 04:38:25 PM CEST 2022
*/
#include <iostream>
#include <string.h>
#include <string>

#include <mpi.h> 

#include "solvers.h"

#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif

#ifdef USE_ASCENT
#include "AscentAdaptor.h"
#endif

int main(int argc, char *argv[])
{
  int grid_resolution=64;
  int Catalyst_argc = argc;
  std::set<std::string> meshtypes = {"uniform", "rectilinear", "structured", "unstructured"};
  std::string meshtype = "uniform";
    
  // we first capture the only two possible arguments allowed --res= and --mesh=
  for (auto cc = 1; cc < argc; ++cc)
  {
    if (strncmp(argv[cc], "--res=", 6) == 0)
    {
      auto n = atoi(&argv[cc][6]);
      if (n > 0 && n < 1026)
        grid_resolution = n;
      else
        grid_resolution = 64;
      Catalyst_argc--;
    }
    else if (strncmp(argv[cc], "--mesh=", 7) == 0)
    {
      auto pos = meshtypes.find(&argv[cc][7]);
      if(pos != meshtypes.end())
        meshtype = &argv[cc][7]; // this is a supported meshtype
      else
        {
        std::cerr  << "mesh type not implemented\n";
        std::cerr << "Supported types: ";
        for (auto &m : meshtypes)
          std::cerr << m << ", ";
        std::cerr << std::endl;
        exit(1);
        }
      Catalyst_argc--;
    }
  }

  std::cout  << "Creating mesh of type "<< meshtype << " of resolution " << grid_resolution << "x" << grid_resolution << std::endl;
  simulation_data sim = {.resolution = grid_resolution, .mesh = meshtype};
  SimInitialize(&sim);

  sim.cart_dims[0] = sim.cart_dims[1] = 0;
  int PartitioningDimension = 2; // want a 2D MPI partitioning. otherwise set to 1.
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &sim.par_rank); /* get current process id */
  MPI_Comm_size(MPI_COMM_WORLD, &sim.par_size); /* get # procs from env or */

  MPI_Partition(PartitioningDimension, &sim);

  neighbors(&sim);

// We use (bx + 2) grid points in the X direction, i.e. interior points plus 2 b.c. points
// We use (by + 2) grid points in the Y direction, i.e. interior points plus 2 b.c. points
// decompose the domain

  AllocateGridMemory(&sim);

  set_initial_bc(&sim);

#ifdef USE_CATALYST
  CatalystAdaptor::Initialize(Catalyst_argc, &argv[argc-Catalyst_argc]);
  std::cout << "CatalystInitialize" << std::endl;
#endif
#ifdef USE_ASCENT
  AscentAdaptor::Initialize(Catalyst_argc, &argv[argc-Catalyst_argc], &sim);
  std::cout << "AscentInitialize" << std::endl;
#endif

  while (sim.gdel > TOL)
    {
    simulate_one_timestep(&sim);

#ifdef USE_CATALYST
    CatalystAdaptor::Execute(sim); // sim.iter*0.1, sim.iter*0.1, grid, attributes);
#endif
#ifdef USE_ASCENT
    AscentAdaptor::Execute(sim); // sim.iter*0.1, sim.iter*0.1, grid, attributes);
#endif
    }

  if (!sim.par_rank)
    std::cout << "Stopped at iteration " << sim.iter << " . Maximum error = " << sim.gdel << std::endl;

  WriteFinalGrid(&sim);

#ifdef USE_CATALYST
  CatalystAdaptor::Finalize();
#endif

#ifdef USE_ASCENT
  AscentAdaptor::Finalize();
#endif

  FreeGridMemory(&sim);

  MPI_Finalize();

  return (0);
}

