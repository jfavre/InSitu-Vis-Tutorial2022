#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include "solvers.h"
#include <catalyst.hpp>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace CatalystAdaptor
{
static std::vector<std::string> filesToValidate;

/**
 * In this example, we show how we can use Catalysts's C++
 * wrapper around conduit's C API to create Conduit nodes.
 * This is not required. A C++ adaptor can just as
 * conveniently use the Conduit C API to setup the
 * `conduit_node`. However, this example shows that one can
 * indeed use Catalyst's C++ API, if the developer so chooses.
 */
void Initialize(int argc, char* argv[])
{
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

      // alternatively, use this form to pass optional parameters to the script.
      const auto name = "catalyst/scripts/script" + std::to_string(cc - 1);
      node[name + "/filename"].set_string(path);
    }
  }

  // indicate that we want to load ParaView-Catalyst
  node["catalyst_load/implementation"].set_string("paraview");
  node["catalyst_load/search_paths/paraview"] = PARAVIEW_IMPL_DIR;

  catalyst_status err = catalyst_initialize(conduit_cpp::c_node(&node));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to initialize Catalyst: " << err << std::endl;
  }
}

void Execute(simulation_data& sim) //int cycle, double time, Grid& grid, Attributes& attribs)
{
  conduit_cpp::Node exec_params;

  // add time/cycle information
  auto state = exec_params["catalyst/state"];
  state["timestep"].set(sim.iter);
  state["time"].set(sim.iter*0.1);

  // Add channels.
  // We only have 1 channel here. Let's name it 'grid'.
  auto channel = exec_params["catalyst/channels/grid"];

  // Since this example is using Conduit Mesh Blueprint to define the mesh,
  // we set the channel's type to "mesh".
  channel["type"].set("mesh");

  // now create the mesh.
  auto mesh = channel["data"];

  if(sim.mesh == "rectilinear")
    {
    mesh["coordsets/coords/values/x"].set_external(sim.cx, (sim.bx + 2));
    mesh["coordsets/coords/values/y"].set_external(sim.cy, (sim.by + 2));
    mesh["coordsets/coords/values/z"].set(0.0);
    mesh["coordsets/coords/type"].set(sim.mesh);
    }
  else if(sim.mesh == "uniform")
    {
    //std::cout << "Uniform Grid dimensions =[" << (sim.local_extents[1] - sim.local_extents[0] + 1) << ", " << (sim.local_extents[3] - sim.local_extents[2] + 1) << ", 1]"<< std::endl;
    mesh["coordsets/coords/dims/i"].set(sim.local_extents[1] - sim.local_extents[0] + 1);
    mesh["coordsets/coords/dims/j"].set(sim.local_extents[3] - sim.local_extents[2] + 1);
    mesh["coordsets/coords/dims/k"].set(1);
    
    //std::cout << "Uniform Grid Origin =[" << sim.cx[0] << ", " << sim.cy[0] << ", 0.]"<< std::endl;
    mesh["coordsets/coords/origin/x"].set(sim.cx[0]);
    mesh["coordsets/coords/origin/y"].set(sim.cy[0]);
    mesh["coordsets/coords/origin/z"].set(0.0);
    mesh["coordsets/coords/type"].set(sim.mesh);

    float spacing = 1.0/(sim.resolution+1.0);
    mesh["coordsets/coords/spacing/dx"].set(spacing);
    mesh["coordsets/coords/spacing/dy"].set(spacing);
    mesh["coordsets/coords/spacing/dz"].set(spacing);
    }

  else if((sim.mesh == "structured") || (sim.mesh == "unstructured"))
    {
    //std::cout << "Explicit Grid dimensions =[" << (sim.bx + 2) * (sim.by + 2) << ", " << (sim.bx + 2) * (sim.by + 2) << ", " << (sim.bx + 2) * (sim.by + 2)<< std::endl;
    mesh["coordsets/coords/type"].set("explicit");
    mesh["coordsets/coords/values/x"].set_external(sim.explicit_cx, (sim.bx + 2) * (sim.by + 2),0,sizeof(float));
    mesh["coordsets/coords/values/y"].set_external(sim.explicit_cy, (sim.bx + 2) * (sim.by + 2),0,sizeof(float));
    mesh["coordsets/coords/values/z"].set_external(sim.explicit_cz, (sim.bx + 2) * (sim.by + 2),0,sizeof(float));
    }

  // add topology.
    mesh["topologies/mesh/type"].set(sim.mesh);
    mesh["topologies/mesh/coordset"].set("coords");

    if(sim.mesh == "unstructured")
      {
      mesh["topologies/mesh/elements/shape"] = "quad";
      mesh["topologies/mesh/elements/connectivity"].set_external(sim.connectivity, 4 * (sim.bx + 1) * (sim.by + 1));
      }
    else if(sim.mesh == "structured")
      {
      mesh["topologies/mesh/elements/dims/i"].set(sim.local_extents[1] - sim.local_extents[0] + 1);
      mesh["topologies/mesh/elements/dims/j"].set(sim.local_extents[3] - sim.local_extents[2] + 1);
      mesh["topologies/mesh/elements/dims/k"].set(1);
      }
      
  // Finally, add fields.
  auto fields = mesh["fields"];
  // temperature is vertex-data.
  fields["temperature/association"].set("vertex");
  fields["temperature/type"].set("scalar");
  fields["temperature/topology"].set("mesh");
  fields["temperature/volume_dependent"].set("false");
  // Conduit supports zero copy, allowing a Conduit Node to describe and
  // point to externally allocated data
  fields["temperature/values"].set_external(sim.Temp, (sim.bx + 2) * (sim.by + 2));
/*
  fields["vtkGhostType/association"].set("vertex");
  fields["vtkGhostType/topology"].set("mesh");
  fields["vtkGhostType/volume_dependent"].set("false");
  fields["vtkGhostType/values"].set_external(sim.Ghost, (sim.bx + 2) * (sim.by + 2));
  */
  catalyst_status err = catalyst_execute(conduit_cpp::c_node(&exec_params));
  if (err != catalyst_status_ok)
  {
    std::cerr << "ERROR: Failed to execute Catalyst: " << err << std::endl;
  }
}

void Finalize()
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
