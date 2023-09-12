/*Module demonstrating the use of Ascent and Catalyst for in-situ visualization"""
##############################################################################
# A simple simulator for a 2D problem, with an in-situ coupling
# The data generation parameters for the vector field
# come from https://shaddenlab.berkeley.edu/uploads/LCS-tutorial/examples.html
# with the Ascent library   https://ascent.readthedocs.io/en/latest/#
# with the Catalyst library https://catalyst-in-situ.readthedocs.io/en/latest/index.html
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this serial version runs until completion and saves images of the scalar field
# at regular intervals. It saves the final solution to a blueprint HDF5 file
#
# Tested Thu 27 Apr 10:36:20 CEST 2023
#
*/
#include <math.h>
#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>

#include "double_gyre.h"
using double_gyre::simulation;
  
#ifdef USE_CATALYST
#include "double_gyre_catalyst.h"
using namespace CatalystAdaptor;
#endif

#ifdef USE_ASCENT
#include "double_gyre_ascent.h"
using namespace AscentAdaptor;
#endif

int main(int argc, char **argv)
{
  if(argc < 4 || argc > 5){
    std::cerr << "Syntax: double_gyre_ascent x-resolution y-resolution nb_timesteps\n";
    exit(1);
  }
  simulation.AllocateGrid(atoi(argv[1]), atoi(argv[2])); // x-resolution and y-resolution
  int max_iterations = atoi(argv[3]);

#ifdef USE_ASCENT
  Ascent_Initialize();
#endif
#ifdef USE_CATALYST
  Catalyst_Initialize(1, &argv[4]);
#endif

  for(auto iteration=0; iteration < max_iterations; iteration++)
    {
    simulation.compute_step();
#ifdef USE_ASCENT
    Ascent_Execute(10); //frequency = 10. Execute 1 every 10 iterations
#endif
#ifdef USE_CATALYST
    Catalyst_Execute();
#endif
    }
  
#ifdef USE_ASCENT
  Ascent_Finalize();
#endif
#ifdef USE_CATALYST
  Catalyst_Finalize();
#endif
  return 0;
}
