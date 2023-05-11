#include <math.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "double_gyre_ascent.h"

#include "double_gyre.h"

#include <ascent/ascent.hpp>
#include "conduit_blueprint.hpp"

namespace AscentAdaptor
{
  ascent::Ascent ascent;
  conduit::Node mesh;
  conduit::Node actions; // default actions can also be overidden by file ascent_actions.yaml

  using double_gyre::simulation;
  
void Ascent_Initialize()
{
  ascent.open();
  mesh["coordsets/coords/type"] = "uniform";
  mesh["coordsets/coords/dims/i"] = simulation.xres;
  mesh["coordsets/coords/dims/j"] = simulation.yres;

  mesh["topologies/mesh/type"] = "uniform";
  mesh["topologies/mesh/coordset"] = "coords";

  mesh["coordsets/coords/origin/x"] = simulation.grid_bounds[0];
  mesh["coordsets/coords/origin/y"] = simulation.grid_bounds[2];
  mesh["coordsets/coords/spacing/dx"] = simulation.grid_bounds[1] / (simulation.xres - 1.0);
  mesh["coordsets/coords/spacing/dy"] = simulation.grid_bounds[3] / (simulation.yres - 1.0);

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

// verify the mesh we created conforms to the blueprint
  conduit::Node verify_info;
  if(!conduit::blueprint::mesh::verify(mesh, verify_info))
    std::cerr << "DoubleGyre Mesh Verify failed!" << std::endl;
  else{
    std::cerr << "DoubleGyre Mesh verify success!" << std::endl;
   //print(self.mesh.to_yaml())
  }
  conduit::Node &add_act0 = actions.append();
  add_act0["action"] = "add_pipelines";
  conduit::Node &pipelines = add_act0["pipelines"];
  pipelines["pl1/f1/type"] = "vector_magnitude";
  pipelines["pl1/f1/params/field"] = "Velocity";
  pipelines["pl1/f1/params/output_name"] = "velocity_mag2d";

  pipelines["pl2/f1/type"] = "vorticity";
  pipelines["pl2/f1/params/field"] = "Velocity";
  pipelines["pl2/f1/params/output_name"] = "Mvorticity";
        
  pipelines["pl2/f2/type"] = "vector_magnitude";
  pipelines["pl2/f2/params/field"] = "Mvorticity";
  pipelines["pl2/f2/params/output_name"] = "vorticity_mag";

  conduit::Node &add_act1 = actions.append();
  add_act1["action"] = "add_scenes";

  conduit::Node &scenes = add_act1["scenes"];
  scenes["s1/plots/p1/type"] = "pseudocolor";
  scenes["s1/plots/p1/pipeline"] = "pl1";
  scenes["s1/plots/p1/field"] = "velocity_mag2d";
  scenes["s1/image_prefix"] = "vel_mag.%04d";

  scenes["s2/plots/p1/type"] = "pseudocolor";
  scenes["s2/plots/p1/pipeline"] = "pl2";
  scenes["s2/plots/p1/field"] = "vorticity_mag";
  scenes["s2/image_prefix"] = "vort_mag.%04d";
  std::cout << actions.to_yaml() << std::endl;
};

void Ascent_Execute(int frequency)
{
  if(!(simulation.iteration % frequency)){
    mesh["state/cycle"] = simulation.iteration;
    ascent.publish(mesh);
    ascent.execute(actions);
  }
};

void Ascent_Finalize()
{
  ascent.close();
};

};
