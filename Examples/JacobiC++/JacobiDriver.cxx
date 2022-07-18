/*
 A parallel Jacobi solver for the Laplacian equation in 2D
 Written by Jean M. Favre, Swiss National Supercomputing Center
 Last tested Thu Oct 28 01:26:47 PM CEST 2021

*/
#include <iostream>
#include <string>
#ifdef PARALLEL
#include <mpi.h> 
#endif

#include "solvers.h"

#ifdef USE_CATALYST
#include "CatalystAdaptor.h"
#endif

int main(int argc, char *argv[])
{
  int grid_resolution;
  std::string meshtype;
  int Catalyst_argc = argc;
  std::set<std::string> meshtypes = {"uniform", "rectilinear", "structured", "unstructured"};
    
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

#ifdef PARALLEL
  sim.cart_dims[0] = sim.cart_dims[1] = 0;
  int PartitioningDimension = 2; // want a 2D MPI partitioning. otherwise set to 1.
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &sim.par_rank); /* get current process id */
  MPI_Comm_size(MPI_COMM_WORLD, &sim.par_size); /* get # procs from env or */

  MPI_Partition(PartitioningDimension, &sim);

  neighbors(&sim);
#endif

// We use (bx + 2) grid points in the X direction, i.e. interior points plus 2 b.c. points
// We use (by + 2) grid points in the Y direction, i.e. interior points plus 2 b.c. points
  // decompose the domain

  AllocateGridMemory(&sim);

  set_initial_bc(&sim);

#ifdef USE_CATALYST
  CatalystAdaptor::Initialize(Catalyst_argc, &argv[argc-Catalyst_argc]);
  std::cout << "CatalystInitialize" << std::endl;
#endif

  while (sim.gdel > TOL)
    {
    simulate_one_timestep(&sim);

#ifdef USE_CATALYST
    CatalystAdaptor::Execute(sim); // sim.iter*0.1, sim.iter*0.1, grid, attributes);
#endif
    }

#ifdef PARALLEL
  if (!sim.par_rank)
#endif
  std::cout << "Stopped at iteration " << sim.iter << " . Maximum error = " << sim.gdel << std::endl;

  WriteFinalGrid(&sim);

#ifdef USE_CATALYST
  CatalystAdaptor::Finalize();
#endif

  FreeGridMemory(&sim);

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return (0);
}

