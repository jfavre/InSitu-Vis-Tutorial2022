#ifndef AscentAdaptor_h
#define AscentAdaptor_h

#include "solvers.h"
#include <ascent/ascent.hpp>
#include "conduit_blueprint.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace AscentAdaptor
{
  ascent::Ascent ascent;
  conduit::Node mesh;
  conduit::Node actions; // default actions can also be overidden by file ascent_actions.yaml

void Initialize(int argc, char* argv[], const simulation_data *sim)
{
  std::cout << "AscentInitialize.........................................\n";
  conduit::Node ascent_options;
  ascent_options["mpi_comm"] = MPI_Comm_c2f(MPI_COMM_WORLD);
  ascent.open(ascent_options);
  
  if(sim->mesh == "rectilinear")
    {
    mesh["coordsets/coords/values/x"].set_external(sim->cx, (sim->bx + 2));
    mesh["coordsets/coords/values/y"].set_external(sim->cy, (sim->by + 2));
    mesh["coordsets/coords/type"].set(sim->mesh);
    }
  else if(sim->mesh == "uniform")
    {
    //std::cout << "Uniform Grid dimensions =[" << (sim.local_extents[1] - sim.local_extents[0] + 1) << ", " << (sim.local_extents[3] - sim.local_extents[2] + 1) << ", 1]"<< std::endl;
    mesh["coordsets/coords/dims/i"].set(sim->local_extents[1] - sim->local_extents[0] + 1);
    mesh["coordsets/coords/dims/j"].set(sim->local_extents[3] - sim->local_extents[2] + 1);
    // do not specify the 3rd dimension with a dim of 1, a z_origin, and a z_spacing
    
    //std::cout << "Uniform Grid Origin =[" << sim.cx[0] << ", " << sim.cy[0] << ", 0.]"<< std::endl;
    mesh["coordsets/coords/origin/x"].set(sim->cx[0]);
    mesh["coordsets/coords/origin/y"].set(sim->cy[0]);
    mesh["coordsets/coords/type"].set(sim->mesh);

    float spacing = 1.0/(sim->resolution+1.0);
    mesh["coordsets/coords/spacing/dx"].set(spacing);
    mesh["coordsets/coords/spacing/dy"].set(spacing);
    }

  else if((sim->mesh == "structured") || (sim->mesh == "unstructured"))
    {
    //std::cout << "Explicit Grid dimensions =[" << (sim.bx + 2) * (sim.by + 2) << ", " << (sim.bx + 2) * (sim.by + 2) << ", " << (sim.bx + 2) * (sim.by + 2)<< std::endl;
    mesh["coordsets/coords/type"].set("explicit");
    mesh["coordsets/coords/values/x"].set_external(sim->explicit_cx, (sim->bx + 2) * (sim->by + 2), 0, sizeof(double));
    mesh["coordsets/coords/values/y"].set_external(sim->explicit_cy, (sim->bx + 2) * (sim->by + 2), 0, sizeof(double));
    }

  // add topology.
    mesh["topologies/mesh/type"].set(sim->mesh);
    mesh["topologies/mesh/coordset"].set("coords");

    if(sim->mesh == "unstructured")
      {
      mesh["topologies/mesh/elements/shape"] = "quad";
      mesh["topologies/mesh/elements/connectivity"].set_external(sim->connectivity, 4 * (sim->bx + 1) * (sim->by + 1));
      }
    else if(sim->mesh == "structured")
      {
      mesh["topologies/mesh/elements/dims/i"].set(sim->local_extents[1] - sim->local_extents[0]);
      mesh["topologies/mesh/elements/dims/j"].set(sim->local_extents[3] - sim->local_extents[2]);
      }
      
  // temperature is vertex-data.
  mesh["fields/temperature/association"].set("vertex");
  mesh["fields/temperature/type"].set("scalar");
  mesh["fields/temperature/topology"].set("mesh");
  mesh["fields/temperature/volume_dependent"].set("false");
  mesh["fields/temperature/values"].set_external(sim->Temp, (sim->bx + 2) * (sim->by + 2));
  
  conduit::Node verify_info;
  if (!conduit::blueprint::mesh::verify(mesh, verify_info))
    {
    // verify failed, print error message
    CONDUIT_INFO("blueprint verify failed!" + verify_info.to_json());
    }
  else CONDUIT_INFO("blueprint verify success!" + verify_info.to_json());

  conduit::Node &add_action = actions.append();
  
  add_action["action"] = "add_scenes";
  conduit::Node &scenes       = add_action["scenes"];
  scenes["view/plots/p1/type"]  = "pseudocolor";
  scenes["view/plots/p1/field"] = "temperature";
  scenes["view/image_prefix"] = "view_%04d";

  //std::cout << mesh.to_yaml();
}

void Execute(simulation_data& sim) //int cycle, double time, Grid& grid, Attributes& attribs)
{
  mesh["state/cycle"].set(sim.iter);
  mesh["state/time"].set(sim.iter*0.1);

  ascent.publish(mesh);
  ascent.execute(actions);
}

void Finalize()
{
  ascent.close();
  std::cout << "AscentFinalize.........................................\n";
}
}
#endif
