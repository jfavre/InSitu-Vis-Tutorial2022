#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.hpp>
#include <catalyst_conduit_blueprint.hpp>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "double_gyre_catalyst.h"

#include "double_gyre.h"

namespace CatalystAdaptor
{
  using double_gyre::simulation;
  static std::vector<std::string> filesToValidate;
  conduit_cpp::Node exec_params;

void Catalyst_Initialize(int argc, char* argv[])
{
  std::cout << "CatalystInitialize" << std::endl;

  conduit_cpp::Node node;
  for (int cc = 0; cc < argc; ++cc)
  {
    if (strcmp(argv[cc], "--output") == 0 && (cc + 1) < argc)
    {
      node["catalyst/pipelines/0/type"].set("io");
      node["catalyst/pipelines/0/filename"].set(argv[cc + 1]);
      node["catalyst/pipelines/0/channel"].set("grid");
      ++cc;
    }
    else if (strcmp(argv[cc], "--exists") == 0 && (cc + 1) < argc)
    {
      filesToValidate.push_back(argv[cc + 1]);
      ++cc;
    }
    else
    {
      const auto path = std::string(argv[cc]);
      // note: one can simply add the script file as follows:
      // node["catalyst/scripts/script" + std::to_string(cc - 1)].set_string(path);
    std::cout << "Using PV Python script : " << path << std::endl;
      // alternatively, use this form to pass optional parameters to the script.
      const auto name = "catalyst/scripts/script" + std::to_string(cc - 1);
      node[name + "/filename"].set_string(path);
    }
  }

  // indicate that we want to load ParaView-Catalyst
  node["catalyst_load/implementation"].set_string("paraview");
  // the env variable CATALYST_IMPLEMENTATION_PATHS should indicate where to find the ParaView specific implementation
  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to initialize Catalyst: " << err << std::endl;
  }
  
  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Since this example is using Conduit Mesh Blueprint to define the mesh,
  // we set the channel's type to "mesh".
  channel["type"].set("mesh");

  // now create the mesh.
  auto mesh = channel["data"];
  mesh["coordsets/coords/type"] = "uniform";
  mesh["coordsets/coords/dims/i"] = simulation.xres;
  mesh["coordsets/coords/dims/j"] = simulation.yres;
  mesh["coordsets/coords/dims/k"] = 1;
  
  mesh["topologies/mesh/type"] = "uniform";
  mesh["topologies/mesh/coordset"] = "coords";

  mesh["coordsets/coords/origin/x"] = simulation.grid_bounds[0];
  mesh["coordsets/coords/origin/y"] = simulation.grid_bounds[2];
  mesh["coordsets/coords/origin/z"] = 0.;
    
  mesh["coordsets/coords/spacing/dx"] = simulation.grid_bounds[1] / (simulation.xres - 1.0);
  mesh["coordsets/coords/spacing/dy"] = simulation.grid_bounds[3] / (simulation.yres - 1.0);
  mesh["coordsets/coords/spacing/dz"] = simulation.grid_bounds[3] / (simulation.yres - 1.0);
 
  mesh["fields/vel_x/association"] = "vertex";
  mesh["fields/vel_x/topology"] = "mesh";
  mesh["fields/vel_x/values"].set_external(simulation.vel_x);

  mesh["fields/vel_y/association"] = "vertex";
  mesh["fields/vel_y/topology"] = "mesh";
  mesh["fields/vel_y/values"].set_external(simulation.vel_y);

  mesh["fields/Velocity/association"] = "vertex";
  mesh["fields/Velocity/topology"] = "mesh";
  mesh["fields/Velocity/values/u"].set_external(simulation.vel_x);
  mesh["fields/Velocity/values/v"].set_external(simulation.vel_y);
  mesh["fields/Velocity/values/w"].set_external(simulation.vel_z);

  conduit_cpp::Node verify_info;
  if(!conduit_cpp::BlueprintMesh::verify(mesh, verify_info))
    {
    std::cerr << "DoubleGyre Mesh Verify failed!" << std::endl;
    std::cerr << verify_info.to_yaml()<< std::endl;
    }
}

void Catalyst_Execute()
{
  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(simulation.iteration);
  state["time"].set(simulation.iteration*0.1);
  
  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to execute Catalyst: " << err << std::endl;
  }
}

void Catalyst_Finalize()
{
  conduit_cpp::Node node;
  catalyst_status err = catalyst_finalize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to finalize Catalyst: " << err << std::endl;
  }

  for (const auto& fname : filesToValidate)
  {
    std::ifstream istrm(fname.c_str(), std::ios::binary);
    if (!istrm.is_open())
    {
      std::cerr << "ERROR: Failed to open file '" << fname.c_str() << "'." << std::endl;
    }
  }
}
}

#endif
