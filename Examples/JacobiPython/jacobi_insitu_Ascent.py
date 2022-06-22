##############################################################################
# A simple simulator for the heat equation in 2D, with an in-situ coupling
# with the Ascent library https://ascent.readthedocs.io/en/latest/#
#
# Author: Jean M. Favre, Swiss National Supercomputing Center
#
# this serial version runs until completion and saves images of the scalar field
# at regular intervals. It saves the final solution to a blueprint HDF5 file
#
# run: python3 jacobi_insitu_Ascent.py
#
# Tested with Python 3.9.7, Wed Apr 27 12:02:12 PM CEST 2022
#
# module load daint-gpu Ascent matplotlib
#
##############################################################################
import sys, math
import numpy as np
import conduit
import conduit.blueprint
import ascent
import matplotlib.pyplot as plt

class Simulation:
    """
    A simple 4-point stencil simulation for the heat equation
    
    Attributes
    ----------
    resolution : int
        the number of grid points on the I and J axis (default 64)
    iterations : int
        the maximum number of iterations (default 100)
    """
    def __init__(self, resolution=64, iterations=100):
        self.iteration = 0 # current iteration
        self.Max_iterations = iterations
        self.xres = resolution
        self.yres = self.xres
        self.dx = 1.0 / (self.xres + 1)

    def Initialize(self):
        """ 2 additional boundary points are added. Iterations will only touch
        the internal grid points.
        """
        self.rmesh_dims = [self.yres + 2, self.xres + 2]
        #print("grid dimensions = ", self.rmesh_dims)
        self.v = np.zeros(self.rmesh_dims)
        self.vnew = np.zeros([self.yres, self.xres])
        self.set_initial_bc()

    def set_initial_bc(self):
        """ initial values set to 0 except on bottom and top wall """
        #first (bottom) row
        self.v[0,:] = [math.sin(math.pi * j * self.dx)
                       for j in range(self.rmesh_dims[1])]
        #last (top) row
        self.v[-1,:] = self.v[0,:]* math.exp(-math.pi)

    def Finalize(self):
        """plot the scalar field iso-contour lines"""
        fig, ax = plt.subplots()
        CS = ax.contour(self.v, levels=10)
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('Temperature iso-contours')
        plt.savefig('Temperature-iso-contours.01.png')
        #plt.show()
        
    def SimulateOneTimestep(self):
        self.iteration += 1
        # print("Simulating time step: iteration=%d" % self.iteration)

        self.vnew = 0.25 * ( self.v[2:, 1:-1]  + # north neighbor
                             self.v[0:-2, 1:-1] + # south neighbor
                             self.v[1:-1, 2:] + # east neighbor
                             self.v[1:-1, :-2]) # west neighbor
        # copy vnew to the interior region of v, leaving the boundary walls untouched.
        self.v[1:-1,1:-1] = self.vnew.copy()

    def MainLoop(self):
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()

# we now define a sub-class of Simulation to add a Conduit node and Ascent action

class Simulation_With_Ascent(Simulation):
    def __init__(self, resolution=64, iterations=100, meshtype="uniform"):
        Simulation.__init__(self, resolution, iterations)
        self.MeshType = meshtype
        if meshtype == "rectilinear":
          self.xc = np.linspace(0, 1, self.xres + 2)
          self.yc = np.linspace(0, 1, self.yres + 2)
        else:
          if meshtype == "structured" or meshtype == "unstructured":
            self.xc, self.yc = np.meshgrid(np.linspace(0, 1, self.xres + 2),
                                           np.linspace(0, 1, self.yres + 2),
                                           indexing='xy')
        if meshtype == "unstructured":
          self.conn = np.zeros(((self.xres + 1) * (self.yres + 1) * 4), dtype=np.int32)
          i=0
          for iy in range(self.yres+1):
            for ix in range(self.xres+1):
              self.conn[4*i+0] = ix + iy*(self.xres + 2)
              self.conn[4*i+1] = ix + (iy+1)*(self.xres + 2)
              self.conn[4*i+2] = ix + (iy+1)*(self.xres + 2)+ 1
              self.conn[4*i+3] = ix + iy*(self.xres + 2) + 1
              i += 1
    # Add Ascent mesh definition
    def Initialize(self):
        Simulation.Initialize(self)
        # set options to allow errors propagate to python
        ascent_opts = conduit.Node()
        ascent_opts["exceptions"] = "forward"
        # open Ascent
        self.a = ascent.Ascent()
        self.a.open(ascent_opts)

        # setup a uniform mesh
        self.mesh = conduit.Node()

        # create the coordinate set
        if self.MeshType == "structured" or self.MeshType == "unstructured":
          self.mesh["coordsets/coords/type"] = "explicit"
        else:
          self.mesh["coordsets/coords/type"] = self.MeshType

        self.mesh["coordsets/coords/dims/i"] = self.xres + 2
        self.mesh["coordsets/coords/dims/j"] = self.yres + 2
        self.mesh["topologies/mesh/type"] = self.MeshType
        self.mesh["topologies/mesh/coordset"] = "coords"

        if self.MeshType == "uniform":
          # add origin and spacing to the coordset (optional)
          self.mesh["coordsets/coords/origin/x"] = 0.0
          self.mesh["coordsets/coords/origin/y"] = 0.0
          self.mesh["coordsets/coords/spacing/dx"] = self.dx
          self.mesh["coordsets/coords/spacing/dy"] = self.dx
        else:
          self.mesh["coordsets/coords/values/x"].set_external(self.xc.ravel())
          self.mesh["coordsets/coords/values/y"].set_external(self.yc.ravel())

        if self.MeshType == "structured":
          self.mesh["topologies/mesh/elements/dims/i"] = np.int32(self.xres + 1)
          self.mesh["topologies/mesh/elements/dims/j"] = np.int32(self.yres + 1)

        if self.MeshType == "unstructured":
          self.mesh["topologies/mesh/elements/shape"] = "quad"
          self.mesh["topologies/mesh/elements/connectivity"].set_external(self.conn)
        # create a vertex associated field called "temperature"
        self.mesh["fields/temperature/association"] = "vertex";
        self.mesh["fields/temperature/topology"] = "mesh";
        # set_external does not handle multidimensional numpy arrays or
        # multidimensional complex strided views into numpy arrays.
        # Views that are effectively 1D-strided are supported.
        self.mesh["fields/temperature/values"].set_external(self.v.ravel())

        # make sure the mesh we created conforms to the blueprint
        verify_info = conduit.Node()
        if not conduit.blueprint.mesh.verify(self. mesh,verify_info):
          print("Mesh Verify failed!")
          print(verify_info.to_yaml())
        else:
          print("Mesh verify success!")
    
        # print the mesh we created
        # print(self.mesh.to_yaml())

        # setup actions
        self.actions = conduit.Node()
        add_act = self.actions.append()
        add_act["action"] = "add_scenes"

        # declare a scene (s1) to render the dataset
        self.scenes = add_act["scenes"]
        # our first scene (named 's1') will render the field 'temperature'
        self.scenes["s1/plots/p1/type"] = "pseudocolor";
        self.scenes["s1/plots/p1/field"] = "temperature";
        # add a second plot to draw the grid lines
        self.scenes["s1/plots/p2/type"] = "mesh";

    def Finalize(self, savedir="./"):
        """ After the final timestep, we save the solution array to disk
        and we close Ascent"""
        Simulation.Finalize(self)
        
        self.a.publish(self.mesh)
        action = conduit.Node()
        add_extr = action.append()
        add_extr["action"] = "add_extracts"
        extracts = add_extr["extracts"]
        extracts["e1/type"]="relay"
        extracts["e1/params/path"] = savedir + "mesh";
        # use HDF5
        extracts["e1/params/protocol"] = "blueprint/mesh/hdf5";
        self.a.execute(action)
        self.a.close()
        
    def MainLoop(self, frequency=100):
      """Ascent actions will be triggered only once every N iterations
    
        Attributes
          ----------
        frequency : int
          the frequency at which to trigger Ascent actions (default 100)
      """
      while self.iteration < self.Max_iterations:
        self.SimulateOneTimestep()
        if not self.iteration % frequency:
          self.mesh["state/cycle"] = self.iteration
          self.scenes["s1/renders/r1/image_name"] = "temperature-ser.%04d" % self.iteration
          # execute the actions
          self.a.publish(self.mesh)
          self.a.execute(self.actions)

def main():
    # sim = Simulation(resolution=64, iterations=500)
    # choices are meshtype="uniform", "rectilinear", "structured", "unstructured"
    sim = Simulation_With_Ascent(resolution=64-2, iterations=500, meshtype="unstructured")
    sim.Initialize()
    sim.MainLoop(frequency=100)
    sim.Finalize(savedir="/mnt/data/")
    
main()

# list all images which have been rendered to disk
import glob
image_files = glob.glob("Temperature*png")
image_files.sort()
print(image_files)


